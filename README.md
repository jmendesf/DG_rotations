# DG_rotations

## Presentation
The program executes a rotation on a 2D image (3D will be implemented later).
The rotation will be interpolated with different methods. 

The goal of this project is to evaluate and analyse each interpolation method.



## Compile and run
cd rotation2D

mkdir build

cd build

cmake ..

make

./rotation2D <path_to_pgm_file>     (still from build directory)

I recommand you use the contourS.pgm file in the sample/directory



## Output
The program will store the output files in output/ file in the rotation2D/ directory (if the program was runned from the build/ directory
