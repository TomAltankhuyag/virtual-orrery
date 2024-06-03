#version 450

#define PI 3.1415926535897932384626433832795
#define NON_INTERSECTING -1.0

#define MAX_TEXTURES  4	
#define SPACE_TEXTURE  0
#define SUN_TEXTURE  1
#define EARTH_TEXTURE  2
#define MOON_TEXTURE  3

// Lighting constants
#define GAMMA 25
#define kd 1
#define ka 0.01 // turn this to 1 to increase brightness
#define ks 1

#define SPACE_RADIUS  20.0 // skybox
#define SUN_RADIUS  0.4
#define EARTH_RADIUS  0.15
#define MOON_RADIUS  0.03

#define TIME_MULTIPLIER 1
#define MOON_DISTANCE_FACTOR 0.4  // distance from earthOrb
#define EARTH_DISTANCE_FACTOR 1.3 // distance from sunOrb

// Bonus implementation constants
#define EARTH_AXIAL_TILT_DEGREES -23.44
#define MOON_AXIAL_TILT_DEGREES 1.54 // in reference to the sun.
#define MOON_ORBITAL_TILT_DEGREES 5.14

layout(location = 0) out vec4 color;

layout(location = 0) in vec3 p;
layout(location = 1) in vec3 d;

layout( push_constant ) uniform constants
{
	mat4 invView; // camera-to-world
	vec4 proj; // (near, far, aspect, fov)
	float time;

} pc;

layout(binding = 0) uniform sampler2D textures[ MAX_TEXTURES ];

vec3 cameraPosition = pc.invView[3].xyz;

struct CelestialOrb {
    vec3 position;
    float radius;
    int textureCode;
    float axialRotation;
    float axialTiltDegree;
};

float TIME = pc.time * TIME_MULTIPLIER;
float LUNAR_PERIOD = TIME / 27.0;
float EARTH_PERIOD = TIME / 365.0;

float earrthX = cos(EARTH_PERIOD) * EARTH_DISTANCE_FACTOR;
float earrthZ = sin(EARTH_PERIOD) * EARTH_DISTANCE_FACTOR;
vec3 earthPosition = vec3(earrthX, 0, earrthZ);

float moonX = (cos(LUNAR_PERIOD) * MOON_DISTANCE_FACTOR) + earthPosition.x;
float moonZ = (sin(LUNAR_PERIOD) * MOON_DISTANCE_FACTOR) + earthPosition.z;
float moonToEarthLength = length(earthPosition - vec3(moonX, 0, moonZ));
vec3 moonPosition = vec3(moonX, moonToEarthLength * tan(radians(MOON_ORBITAL_TILT_DEGREES)), moonZ);

CelestialOrb spaceOrb = {vec3(0), SPACE_RADIUS, SPACE_TEXTURE, 0, 0};
CelestialOrb sunOrb = {vec3(0), SUN_RADIUS, SUN_TEXTURE, LUNAR_PERIOD, 0};
CelestialOrb earthOrb = {earthPosition, EARTH_RADIUS, EARTH_TEXTURE, TIME, EARTH_AXIAL_TILT_DEGREES};
CelestialOrb moonOrb = {moonPosition, MOON_RADIUS, MOON_TEXTURE, LUNAR_PERIOD, MOON_AXIAL_TILT_DEGREES};

CelestialOrb celestialBodies[4] = {
    spaceOrb,
    sunOrb,
    earthOrb,
    moonOrb
};
// Based on the distance from the camera, sort the orbs in descending order
// EXCEPTION: spaceOrb (orbArray[0]) is always rendered first, so it won't be sorted
void descendSort(inout CelestialOrb orbArray[4]) {
    for (int i = 0; i < orbArray.length() - 1; i++) {
        for (int j = 1; j < orbArray.length() - i - 1; j++) {
            if (length(orbArray[j].position - cameraPosition) < length(orbArray[j + 1].position - cameraPosition)) {
                CelestialOrb temp = orbArray[j];
                orbArray[j] = orbArray[j + 1];
                orbArray[j + 1] = temp;
            }
        }
    }
}

float isIntersecting(vec3 rayPosition, vec3 rayVector, float radius ){
    float rayPositionNorm = length(rayPosition);
    float quadraticB = dot(2.0 *rayVector, rayPosition); 
    float discriminant = pow(quadraticB, 2) -4.0*(pow(rayPositionNorm, 2.0) - pow(radius, 2.0));
    if ( discriminant < 0) {
        return NON_INTERSECTING;
    }
    float t1 = 0.5 * (-quadraticB - sqrt(discriminant));
    float t2 = 0.5 * (-quadraticB + sqrt(discriminant));

    float tmin = min(t1, t2);
    float tmax = max(t1, t2);

    if (tmin > 0) {
        return tmin;
    }
    else if (tmax > 0) {
        return tmax;
    }
    else {
        return NON_INTERSECTING;
    }
}
float isShadowIntersecting(CelestialOrb orb, vec3 intersectionPoint) {
    if (orb.textureCode == SUN_TEXTURE) {
        return NON_INTERSECTING;
    }
    CelestialOrb obstacleOrb = orb.textureCode == EARTH_TEXTURE ? moonOrb : earthOrb;

    vec3 lightVector = normalize(vec3(0) - intersectionPoint);
    vec3 surfaceNormal = normalize(intersectionPoint - orb.position);
    float lightAndNormalFactor = max(dot(surfaceNormal, lightVector), 0);
    lightAndNormalFactor = clamp(lightAndNormalFactor, 0, 1);

    if (lightAndNormalFactor > 0) {
        return isIntersecting(intersectionPoint - obstacleOrb.position, lightVector ,obstacleOrb.radius);
    } else {
        return NON_INTERSECTING;
    }
}
vec4 shading(CelestialOrb orb, vec3 intersectionPoint, vec3 textureColor) { 
    if (orb.position == vec3(0)) {
        return vec4(textureColor, 1);
    }
    float t = isShadowIntersecting(orb, intersectionPoint);
    if (t != NON_INTERSECTING) {
        return vec4(textureColor * ka, 1);
    }

    vec3 lightVector = normalize(sunOrb.position - intersectionPoint);
    // Diffuse light
    vec3 surfaceNormal = normalize(intersectionPoint - orb.position);
    float lightAndNormalFactor = max(dot(surfaceNormal, normalize(lightVector)), 0);
    lightAndNormalFactor = clamp(lightAndNormalFactor, 0, 1);
    vec3 diffuseLight = textureColor * lightAndNormalFactor;

    // Specular light
    vec3 viewVector = normalize(cameraPosition - intersectionPoint);
    vec3 specularAngleVector = normalize(lightVector + viewVector);
    float specularAndNormalFactor = dot(surfaceNormal, specularAngleVector);
    specularAndNormalFactor = clamp(specularAndNormalFactor, 0, 1);
    specularAndNormalFactor = pow(specularAndNormalFactor, GAMMA);

    vec3 specularLight = textureColor * specularAndNormalFactor;

    return vec4(((ks * specularLight) + (kd * diffuseLight) + (ka * textureColor)) * textureColor, 1);
}
vec3 positionYAxisUp(vec3 point) {
    point.z = -point.z;
    point = point.xzy;
    return point;
}
vec3 axialTilt(CelestialOrb orb, vec3 point) {
    if (orb.axialTiltDegree == 0) {
       return positionYAxisUp(point);
    }
    float ANGLE = radians(orb.axialTiltDegree);
    float x = point.x * cos(ANGLE) + point.y * sin(ANGLE);
    float y = -point.x * sin(ANGLE) + point.y * cos(ANGLE);

    return positionYAxisUp(vec3(x, y, point.z));
}
void renderOrb(CelestialOrb orb) {
    vec3 rayPosition = p - orb.position;
    vec3 rayVector = normalize(d);

    float t = isIntersecting(rayPosition, rayVector, orb.radius);
    if (t == NON_INTERSECTING) {
        return;
    }
    vec3 intersectionPoint = rayPosition + t*(rayVector);
    vec3 normal = axialTilt(orb, intersectionPoint);

    float phi = acos(normal.z / orb.radius);
    float theta;
    if (abs(normal.x) < 0.001) {
        theta = sign(normal.y)*PI*0.5; 
    } else {
        theta = atan(normal.y, normal.x); 
    }
    vec4 textureColor =  texture(textures[orb.textureCode], vec2(1.0+0.5*(theta + orb.axialRotation)/PI, phi/PI));

    color = shading(orb, p + (t * rayVector), textureColor.xyz);
}

void main() {
    color = vec4(0.4,0,0.4,1.0);
    // Render farthest orb from camera first
    descendSort(celestialBodies);
    for (int i = 0; i < 4; i++) {
        renderOrb(celestialBodies[i]);
    }
}