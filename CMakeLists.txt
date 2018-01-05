cmake_minimum_required(VERSION 2.8)

project(libnmc)

add_library(libnmc STATIC
  src/Complex.cpp
  src/Double.cpp
)
