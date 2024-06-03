/*
 * Copyright 2023 University of Calgary, Visualization and Graphics Group
 */

#include "Object.h"
#include <VulkanLaunchpad.h>
#include <vulkan/vulkan.hpp>
#include <glm/gtx/normal.hpp>
#include "Camera.h"

using std::cout;
using std::endl;
using std::vector;

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

extern vk::PhysicalDevice mPhysicalDevice;
extern vk::Device mDevice;
extern vk::UniqueCommandPool mCommandPool;
extern vk::Queue mQueue;
extern std::string basePath;

// The following variables are related to texture setup
// --- Begin texture setup variables
vector<TextureInfo> textures;
vector<VkSampler> samplers;		// MUST be unique!


VkDescriptorSetLayout descriptorSetLayout;
VkDescriptorPool descriptorPool;
VkDescriptorSet descriptorSet;
// --- End texture setup variables


// buffers that will live on the GPU.
VkBuffer mObjectVertexData;

// A pipeline that can be used for HW4
VkPipeline pipeline;

// Struct to pack object vertex data
// Only positions for the HW4 Starter
struct Vertex {
	glm::vec3 position;
}; 

// Send view and projection information as push constants,
//  as well as any other parameters we need, via this struct
struct ObjectPushConstants {
    glm::mat4 invView;
    glm::vec4 proj;	// {near, far, aspect, fov} rather than the matrix for ray-tracing
    float time;		// you'll use this to place everything in the scene
} pushConstants;

// camera related variables
extern VklCamera* camera;
extern float aspect;
extern float time_since_epoch;
const float near = 0.01f;
const float far = 1000.0f;
const float fov = 45.0f;

// we'll use this in multiple places
vector<std::string> imageFiles = {
	"textures/background.jpg",
	"textures/sun.jpg",
	"textures/earth.jpg",
	"textures/moon.jpg"
	};


// Setup geometry and textures for ray tracing.
void objectCreateGeometryAndBuffers( GLFWwindow* window ) 
{
	// Two triangles that will cover the full screen
	vector<Vertex> vData(6);
	vData[0].position = glm::vec3{ 1.0f, 1.0f, 0.5f};
	vData[1].position = glm::vec3{ 1.0f,-1.0f, 0.5f};
	vData[2].position = glm::vec3{-1.0f, 1.0f, 0.5f};
	vData[3].position = glm::vec3{ 1.0f,-1.0f, 0.5f};
	vData[4].position = glm::vec3{-1.0f,-1.0f, 0.5f};
	vData[5].position = glm::vec3{-1.0f, 1.0f, 0.5f};
	
	const auto device = vklGetDevice();

	cout << "Sending packed vertex-normal data to GPU..." << endl;

	// 1. Vertex BUFFER (Buffer, Memory, Bind 'em together, copy data into it)
	mObjectVertexData = vklCreateHostCoherentBufferAndUploadData(
		vData.data(), sizeof(vData[0]) * vData.size(),
		VK_BUFFER_USAGE_VERTEX_BUFFER_BIT);

	// Now Create the camera and pipeline
	cout << "Now Creating Camera and Pipeline " << endl;
	
	// All texture related setup

	objectCreateTextureImages();
	createTextureImageViews();
	createTextureSamplers();
	createDescriptorPool();
	createDescriptorSet();

	objectCreateCamera( window );
	objectCreatePipeline();

}

// Cleanup buffers and pipeline created on the GPU 
void objectDestroyBuffers() {

	auto device = vklGetDevice();
	vkDeviceWaitIdle( device );
	vklDestroyGraphicsPipeline(pipeline);
	vklDestroyHostCoherentBufferAndItsBackingMemory( mObjectVertexData );

	// cleanup textures
	for( const auto& texture : textures ) {

		vkDestroyImage(device, texture.image, nullptr);
		vkFreeMemory(device, texture.memory, nullptr);
		vkDestroyImageView(device, texture.view, nullptr);
	}
	for( const auto& sampler : samplers )
		vkDestroySampler(device, sampler, nullptr);

	vkDestroyDescriptorPool(device, descriptorPool, nullptr);
	vkDestroyDescriptorSetLayout(device, descriptorSetLayout, nullptr);

	// Cleanup Camera
	vklDestroyCamera( camera);
}

void objectDraw() {
	objectDraw( pipeline );
}

void objectDraw(VkPipeline pipeline)
{
	if (!vklFrameworkInitialized()) {
		VKL_EXIT_WITH_ERROR("Framework not initialized. Ensure to invoke vklFrameworkInitialized beforehand!");
	}
	const vk::CommandBuffer& cb = vklGetCurrentCommandBuffer();
	auto currentSwapChainImageIndex = vklGetCurrentSwapChainImageIndex();
	assert(currentSwapChainImageIndex < vklGetNumFramebuffers());
	assert(currentSwapChainImageIndex < vklGetNumClearValues());

	// use vklCmdBindPipeline for shader hot reloading
	vklCmdBindPipeline(cb, VK_PIPELINE_BIND_POINT_GRAPHICS, pipeline );
	vklBindDescriptorSetToPipeline(descriptorSet, pipeline );
	cb.bindVertexBuffers(0u, { vk::Buffer{ objectGetVertexBuffer() } }, { vk::DeviceSize{ 0 } });
	
    	// upload the matrix to the GPU via push constants
	objectUpdateConstants();
	vklSetPushConstants(
			pipeline, 
			VK_SHADER_STAGE_VERTEX_BIT | VK_SHADER_STAGE_FRAGMENT_BIT, 
			&pushConstants, 
			sizeof(ObjectPushConstants)
		);

	cb.draw(6u, 1u, 0u, 0u);	// it's a full-screen quad, there's only six vertices by definition
}

VkBuffer objectGetVertexBuffer() { return static_cast<VkBuffer>(mObjectVertexData); }

void objectCreatePipeline() {

	// ------------------------------
	// Pipeline creation
	// ------------------------------

	auto const vertShaderName = basePath + std::string("/shaders/starter.vert");
	if( !check_path( vertShaderName ) )
		throw std::runtime_error( "ERROR: Could not load vertex shader!" );
	auto const fragShaderName = basePath + std::string("/shaders/starter.frag");
	if( !check_path( fragShaderName ) )
		throw std::runtime_error( "ERROR: Could not load fragment shader!" );

	VklGraphicsPipelineConfig config{};
	config.enableAlphaBlending = false;

	config.vertexShaderPath = vertShaderName.c_str();
	config.fragmentShaderPath = fragShaderName.c_str();
		
	config.polygonDrawMode = VK_POLYGON_MODE_FILL;
	config.triangleCullingMode = VK_CULL_MODE_NONE;

	// Binding for vertex buffer, using 1 buffer with per-vertex rate.
	// This will send per-vertex data to the GPU.
	config.vertexInputBuffers.emplace_back(VkVertexInputBindingDescription{
		.binding = 0,
		.stride = sizeof(Vertex),
		.inputRate = VK_VERTEX_INPUT_RATE_VERTEX,
	});

	// Positions at locaion 0
	config.inputAttributeDescriptions.emplace_back(VkVertexInputAttributeDescription{
			.location = 0,
			.binding = 0,
			.format = VK_FORMAT_R32G32B32_SFLOAT,
			.offset = offsetof(Vertex, position),
		});

	// all the images we're loading
	config.descriptorLayout.emplace_back( VkDescriptorSetLayoutBinding{
		.binding = 0,
		.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
		.descriptorCount = static_cast<uint32_t>(textures.size()),
		.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT,
		.pImmutableSamplers = nullptr,
	});
            
	// Push constants should be available in both the vertex and fragment shaders
	config.pushConstantRanges.emplace_back(VkPushConstantRange{
			.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT | VK_SHADER_STAGE_VERTEX_BIT,
			.offset = 0,
			.size = sizeof(ObjectPushConstants),
		});
	pipeline = vklCreateGraphicsPipeline( config );
}

// Function to update push constants.
// For the starter example, only the model matrix is updated.
void objectUpdateConstants() {

	pushConstants.invView = glm::inverse(vklGetCameraViewMatrix(camera));	
	pushConstants.proj = glm::vec4{near, far, aspect, fov};
	pushConstants.time = time_since_epoch;
}

// Function to create camera
void objectCreateCamera( GLFWwindow* window ) {
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	// assume a default right handed camera looking down -z.
	
	aspect = static_cast<float>(width) / static_cast<float>(height);

	camera = vklCreateCamera( window, 
		vklCreatePerspectiveProjectionMatrix(glm::radians(fov), 
		aspect, near, far) );
	camera->mPosition = glm::vec3(0.0f, 0.0f, 3.0f);
}

// The code that follows is based on the texture setup code 
// provided by:
// https://vulkan-tutorial.com/Texture_mapping/Images

TextureInfo objectCreateTextureImage( const char* file ) {

	TextureInfo retVal;

	int texWidth, texHeight, texChannels;
	stbi_uc* pixels = stbi_load(file, &texWidth, &texHeight, &texChannels, STBI_rgb_alpha);
	VkDeviceSize imageSize = texWidth * texHeight * STBI_rgb_alpha * sizeof(stbi_uc);

	if (!pixels)
        	throw std::runtime_error("failed to load texture image!");

	VkBuffer stagingBuffer = vklCreateHostCoherentBufferAndUploadData( 
			pixels, 
			static_cast<size_t>(imageSize), 
	 		VK_BUFFER_USAGE_TRANSFER_SRC_BIT );
	
	stbi_image_free(pixels);

	createImage(texWidth, texHeight, VK_FORMAT_R8G8B8A8_SRGB, 
		VK_IMAGE_TILING_OPTIMAL, 
		VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_SAMPLED_BIT,
		VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, 
		retVal.image, retVal.memory);

	transitionImageLayout(retVal.image, VK_FORMAT_R8G8B8A8_SRGB, 
			VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL);
	copyBufferToImage(stagingBuffer, retVal.image, static_cast<uint32_t>(texWidth), static_cast<uint32_t>(texHeight));	
	transitionImageLayout(retVal.image, VK_FORMAT_R8G8B8A8_SRGB, 
			VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL);

	vklDestroyHostCoherentBufferAndItsBackingMemory(stagingBuffer);

	return retVal;
}

void objectCreateTextureImages() {

	// load each image up in turn
	for( const auto& filename : imageFiles ) {
		auto const imagePath = basePath + std::string("/") + filename;
		if( !check_path( imagePath ) )
		        throw std::runtime_error( "ERROR: Could not load texture!" );
		textures.push_back( objectCreateTextureImage(imagePath.c_str()) );
	}
}

void createImage( uint32_t width, uint32_t height, VkFormat format, VkImageTiling tiling, VkImageUsageFlags usage, 
		VkMemoryPropertyFlags properties, VkImage& image, VkDeviceMemory& imageMemory ) {

    VkImageCreateInfo imageInfo{};
    imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;

    imageInfo.imageType = VK_IMAGE_TYPE_2D;
    imageInfo.extent.width = width;
    imageInfo.extent.height = height;
    imageInfo.extent.depth = 1;
    imageInfo.mipLevels = 1;
    imageInfo.arrayLayers = 1;
    imageInfo.format = format;
    imageInfo.tiling = tiling;
    imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    imageInfo.usage = usage;
    imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
    imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    if (vkCreateImage(mDevice, &imageInfo, nullptr, &image) != VK_SUCCESS) {
        throw std::runtime_error("failed to create image!");
    }

    VkMemoryRequirements memRequirements;
    vkGetImageMemoryRequirements(mDevice, image, &memRequirements);

    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;

    allocInfo.allocationSize = memRequirements.size;
    allocInfo.memoryTypeIndex = findMemoryType(memRequirements.memoryTypeBits, properties);

    if (vkAllocateMemory(mDevice, &allocInfo, nullptr, &imageMemory) != VK_SUCCESS) {
        throw std::runtime_error("failed to allocate image memory!");
    }

    vkBindImageMemory(mDevice, image, imageMemory, 0);
}


uint32_t findMemoryType(uint32_t typeFilter, VkMemoryPropertyFlags properties) {
        VkPhysicalDeviceMemoryProperties memProperties;
	vkGetPhysicalDeviceMemoryProperties(mPhysicalDevice, &memProperties);

        for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++) {
            if ((typeFilter & (1 << i)) && (memProperties.memoryTypes[i].propertyFlags & properties) == properties) {
                return i;
            }
        }

        throw std::runtime_error("failed to find suitable memory type!");
}

VkCommandBuffer beginSingleTimeCommands() {

    VkCommandBufferAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;

    allocInfo.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocInfo.commandPool = mCommandPool.get();
    allocInfo.commandBufferCount = 1;

    VkCommandBuffer commandBuffer;
    vkAllocateCommandBuffers(mDevice, &allocInfo, &commandBuffer);

    VkCommandBufferBeginInfo beginInfo{};
    beginInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    beginInfo.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;

    vkBeginCommandBuffer(commandBuffer, &beginInfo);

    return commandBuffer;
}

void endSingleTimeCommands(VkCommandBuffer commandBuffer) {
    vkEndCommandBuffer(commandBuffer);

    VkSubmitInfo submitInfo{};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;

    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;

    vkQueueSubmit(mQueue, 1, &submitInfo, VK_NULL_HANDLE);
    vkQueueWaitIdle(mQueue);

    vkFreeCommandBuffers(mDevice, mCommandPool.get(), 1, &commandBuffer);
}

void copyBuffer(VkBuffer srcBuffer, VkBuffer dstBuffer, VkDeviceSize size) {
    VkCommandBuffer commandBuffer = beginSingleTimeCommands();

    VkBufferCopy copyRegion{};
    copyRegion.size = size;
    vkCmdCopyBuffer(commandBuffer, srcBuffer, dstBuffer, 1, &copyRegion);

    endSingleTimeCommands(commandBuffer);
}

void transitionImageLayout(VkImage image, VkFormat format, VkImageLayout oldLayout, VkImageLayout newLayout) {

	VkCommandBuffer commandBuffer = beginSingleTimeCommands();

	VkImageMemoryBarrier barrier{};
	barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;

	barrier.oldLayout = oldLayout;
	barrier.newLayout = newLayout;

	barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
	barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;

	barrier.image = image;
	barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	barrier.subresourceRange.baseMipLevel = 0;
	barrier.subresourceRange.levelCount = 1;
	barrier.subresourceRange.baseArrayLayer = 0;
	barrier.subresourceRange.layerCount = 1;

	VkPipelineStageFlags sourceStage;
	VkPipelineStageFlags destinationStage;

	if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL) {
    		barrier.srcAccessMask = 0;
    		barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;

    		sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
    		destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
	} else if (oldLayout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL) {
    		barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
    		barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

    		sourceStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
    		destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
	} else {
    		throw std::invalid_argument("unsupported layout transition!");
	}

	vkCmdPipelineBarrier(
    		commandBuffer,
    		sourceStage, destinationStage,
    		0,
    		0, nullptr,
    		0, nullptr,
    		1, &barrier
	);
    endSingleTimeCommands(commandBuffer);
}

void copyBufferToImage(VkBuffer buffer, VkImage image, uint32_t width, uint32_t height) {

   	VkCommandBuffer commandBuffer = beginSingleTimeCommands();

	VkBufferImageCopy region{};
	region.bufferOffset = 0;
	region.bufferRowLength = 0;
	region.bufferImageHeight = 0;

	region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	region.imageSubresource.mipLevel = 0;
	region.imageSubresource.baseArrayLayer = 0;
	region.imageSubresource.layerCount = 1;

	region.imageOffset = {0, 0, 0};
	region.imageExtent = {width, height, 1};

	vkCmdCopyBufferToImage(
    		commandBuffer,
    		buffer,
    		image,
    		VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
    		1,
    		&region
	);

    endSingleTimeCommands(commandBuffer);
}

void createTextureImageViews() {

	// all images use the same basic view, so recycle this struct
	VkImageViewCreateInfo viewInfo{};
	viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;

	viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
	viewInfo.format = VK_FORMAT_R8G8B8A8_SRGB;
	viewInfo.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
	viewInfo.subresourceRange.baseMipLevel = 0;
	viewInfo.subresourceRange.levelCount = 1;
	viewInfo.subresourceRange.baseArrayLayer = 0;
	viewInfo.subresourceRange.layerCount = 1;

	for( auto& texture : textures ) {

		viewInfo.image = texture.image;		// only this values changes
		if (vkCreateImageView(mDevice, &viewInfo, nullptr, &(texture.view)) != VK_SUCCESS) {
	    		throw std::runtime_error("failed to create texture image view!");
		}

	}

}

void createTextureSamplers() {

	// create any samplers we need here
	VkSamplerCreateInfo samplerInfo{};
	samplerInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;

	samplerInfo.magFilter = VK_FILTER_LINEAR;
	samplerInfo.minFilter = VK_FILTER_LINEAR;
	samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_REPEAT;
	samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_REPEAT;
	samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_REPEAT;
	samplerInfo.anisotropyEnable = VK_FALSE; // disabled since we have not talked about anisotropic filtering
	samplerInfo.borderColor = VK_BORDER_COLOR_INT_OPAQUE_BLACK;
	samplerInfo.unnormalizedCoordinates = VK_FALSE; // Use normalized texture coordinates (in the [0,1] range)

	samplerInfo.compareEnable = VK_FALSE;
	samplerInfo.compareOp = VK_COMPARE_OP_ALWAYS;

	// We'll not worry about mipmapping
	samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR;
	samplerInfo.mipLodBias = 0.0f;
	samplerInfo.minLod = 0.0f;
	samplerInfo.maxLod = 0.0f;

	VkSampler sampler;
	if (vkCreateSampler(mDevice, &samplerInfo, nullptr, &sampler) != VK_SUCCESS) {
	       	throw std::runtime_error("failed to create texture sampler!");
	}
	samplers.push_back( sampler );
}

void createDescriptorPool() {

    // tack them onto our pool size pool
    VkDescriptorPoolSize poolTemplate{};
    poolTemplate.type = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
    poolTemplate.descriptorCount = 1u;
    vector<VkDescriptorPoolSize> poolSizes( textures.size(), poolTemplate );

    VkDescriptorPoolCreateInfo poolInfo{};
    poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;

    poolInfo.poolSizeCount = static_cast<uint32_t>(poolSizes.size());
    poolInfo.pPoolSizes = poolSizes.data();
    poolInfo.maxSets = 1u;	// only one frame in flight = only one set needed

    if (vkCreateDescriptorPool(mDevice, &poolInfo, nullptr, &descriptorPool) != VK_SUCCESS) {
        throw std::runtime_error("failed to create descriptor pool!");
    }
}

void createDescriptorSet() {

	// we only need one, to use texture arrays
	VkDescriptorSetLayoutBinding binding = VkDescriptorSetLayoutBinding{
		.binding = 0,
		.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER,
		.descriptorCount = static_cast<uint32_t>( textures.size() ),
		.stageFlags = VK_SHADER_STAGE_FRAGMENT_BIT,
		.pImmutableSamplers = nullptr,
	};
            
	VkDescriptorSetLayoutCreateInfo layoutInfo{};
	layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;

	layoutInfo.bindingCount = 1u;
	layoutInfo.pBindings = &binding;

    	if (vkCreateDescriptorSetLayout(mDevice, &layoutInfo, nullptr, &descriptorSetLayout) != VK_SUCCESS) {
        	throw std::runtime_error("failed to create descriptor set layout!");
    	}

    	VkDescriptorSetAllocateInfo allocInfo{};
	allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;

	allocInfo.descriptorPool = descriptorPool;
	allocInfo.descriptorSetCount = 1u;
	allocInfo.pSetLayouts = &descriptorSetLayout;

	if (vkAllocateDescriptorSets(mDevice, &allocInfo, &descriptorSet) != VK_SUCCESS) {
        	throw std::runtime_error("failed to allocate descriptor sets!");
    	}

	// update all the descriptor sets
	vector<VkDescriptorImageInfo> imageInfo;
	for( uint32_t i = 0; i < textures.size(); i++ ) 

		imageInfo.emplace_back( VkDescriptorImageInfo{
			.sampler = samplers[0],
			.imageView = textures[i].view,
			.imageLayout = VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
			});

	VkWriteDescriptorSet dwrite{};
	dwrite.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;

	dwrite.dstSet = descriptorSet;
	dwrite.dstBinding = 0;
	dwrite.dstArrayElement = 0;
	dwrite.descriptorType = VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER;
	dwrite.descriptorCount = imageInfo.size();
	dwrite.pImageInfo = imageInfo.data();

	vkUpdateDescriptorSets(mDevice, 1u, &dwrite, 0, nullptr);

}

// a helper to check if a file exists and is non-empty
bool check_path( std::string pathName ) {

	const auto path = std::filesystem::path( pathName );
	if( !std::filesystem::exists(path) ) {
		cout << "ERROR: The given file " << path << " does not exist!" <<
			" Please check you gave the path correctly." << endl;
		return false;
		}
	if( std::filesystem::is_empty(path) ) {
		cout << "ERROR: The given file " << path << " is empty!" <<
			" Please check you're pointing to a valid file." << endl;
		return false;
		}

	return true;
	}

