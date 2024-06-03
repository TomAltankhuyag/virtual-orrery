/*
 * Copyright 2023 University of Calgary, Visualization and Graphics Grpup
 */

// stack these to convert a preproccessor string to an actual string
#define L( string ) 		#string
#define TO_LITERAL( string )	L(string)

// if CMake hasn't given us a default assets directory, just guess
#if !defined(ASSET_DIR)
	#define ASSET_DIR ../assets
#endif	

#include <vulkan/vulkan.h>
#include "VulkanLaunchpad.h"
#include <filesystem>
#include <string>

struct TextureInfo {	// a quick helper to track key variables on a per-texture basis

	VkImage image;
	VkDeviceMemory memory;
	VkImageView view;
	VkSampler sampler;		// can be a duplicate of another value

};

void objectCreateGeometryAndBuffers( GLFWwindow* window );
void objectDestroyBuffers();
void objectDraw();

void objectDraw(VkPipeline pipeline);

VkBuffer objectGetVertexBuffer();

void objectCreateTextureImages();
void objectCreateCamera( GLFWwindow* window );
void objectCreatePipeline();
void objectUpdateConstants();

// texture management functions

TextureInfo objectCreateTextureImage( const char* file );
void createImage(uint32_t width, uint32_t height, VkFormat format, VkImageTiling tiling, VkImageUsageFlags usage, VkMemoryPropertyFlags properties, VkImage& image, VkDeviceMemory& imageMemory);
uint32_t findMemoryType(uint32_t typeFilter, VkMemoryPropertyFlags properties);
VkCommandBuffer beginSingleTimeCommands();
void endSingleTimeCommands(VkCommandBuffer commandBuffer);
void copyBuffer(VkBuffer srcBuffer, VkBuffer dstBuffer, VkDeviceSize size);
void transitionImageLayout(VkImage image, VkFormat format, VkImageLayout oldLayout, VkImageLayout newLayout);
void copyBufferToImage(VkBuffer buffer, VkImage image, uint32_t width, uint32_t height);

void createTextureImageViews();
void createTextureSamplers();
void createDescriptorPool();
void createDescriptorSet();

bool check_path( std::string path );
