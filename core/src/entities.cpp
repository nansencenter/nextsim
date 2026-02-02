#include <entities.hpp>

namespace Nextsim
{
namespace entities
{

int getNumVerticesForElementType(int type)
{
    switch(type) {
        case 1:  return 2;  // Line 2
        case 2:  return 3;  // Triangle 3
        case 3:  return 4;  // Quadrangle 4
        case 4:  return 4;  // Tetrahedron 4
        case 5:  return 8;  // Hexahedron 8
        case 8:  return 3;  // Line 3
        case 9:  return 6;  // Triangle 6
        case 15: return 1;  // Point
        case 16: return 8;  // Quadrangle 8
        case 20: return 9;  // Triangle 9
        case 21: return 10; // Triangle 10
        default: return 0;
    }
}

const char* getElementTypeName(int type)
{
    switch(type) {
        case 1:  return "Line 2";
        case 2:  return "Triangle 3";
        case 3:  return "Quadrangle 4";
        case 4:  return "Tetrahedron 4";
        case 5:  return "Hexahedron 8";
        case 8:  return "Line 3";
        case 9:  return "Triangle 6";
        case 15: return "Point";
        case 16: return "Quadrangle 8";
        case 20: return "Triangle 9";
        case 21: return "Triangle 10";
        default: return "Unknown";
    }
}

} // entities
} // Nextsim
