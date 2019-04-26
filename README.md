# DG_rotations
### Presentation
The program executes a rotation of a give angle on a 2D image (3D will be implemented later).
The rotation will be interpolated with different methods. At the moment, nearest neighbor, bilinear and bicubic interpolations are available.

The goal of this project is to evaluate and analyse the results of each interpolation method on various digital shapes, with various topological and geometrical properties.

### Compile and run
From the root directory of the program: 
```
cd rotation2D
mkdir build
cd build
cmake ..
make
./rotation2D <path_to_pgm_file> <angle_in_radians> <_method> [<_lowerThresh><_upperThresh>] 
```

Lower and upper thresholding values are necessary in the case of a non binary image. In the case of a binary image, specifying these two values is not necessary.

The input image must be of .pgm format. 



##### Example :

 ```./rotation2D ../samples/contourS.pgm 1. all 1 135```

``` ./rotation2D ../samples/oval.pgm 1. all``` 
### Output
The program will store the output files in `output/` file in the `rotation2D/` directory (if the program was runned from the `build/` directory. 

The `pre_processing/` folder allows to visualise the steps necessary to the preparation of the data for interpolation.
The `rotation_NN/`, `rotation_BIL/` and `rotation_BIC/` folders contain respectively the results of the rotation using nearest neighbor interpolation, bilinear interpolation and bicubic interpolation, in the .pgm and .eps formats.
