#include "../include/tools.h" 
#include <iostream>

int clampInt(int value, int low, int high)
{
    return value < low ? low : value > high ? high : value; 
}

