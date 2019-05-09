
# DG_rotations
### Presentation
The program executes a rotation of a give angle on a 2D image (3D will be implemented later).
The rotation will be interpolated with different methods. At the moment, nearest neighbor, bilinear and bicubic interpolations are available.

The goal of this project is to evaluate and analyse the results of each interpolation method on various digital shapes, with various topological and geometrical properties.

### Algorithm
The algorithm relies on 3 important steps.

##### Step 1: Distance Map
The input image is binarized if necessary. A distance map is computed for the foreground and the background of the binary image. To do so, the image is duplicated. The duplicate is inversed and a distance map is applied to the duplicate and the original. 
The distance map values of the foreground are set to their opposite value. The two images are then merged together.

##### Step 2: Backward Rotation and Interpolation
An empty image is created with boundaries corresponding to the positions of the 4 corners, forwardly rotated.
The backward rotation is then applied to each point of the destination image. The value of a given point is determined by the interpolation method used. 

 - [Nearest Neighbor](https://en.wikipedia.org/wiki/Nearest-neighbor_interpolation): the rounding function is used. The point takes the value of the corresponding point.
 - [Bilinear Interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation): the four points surrounding the result of the rotation are interpolated and a value is returned.
 - [Bicubic Interpolation](https://en.wikipedia.org/wiki/Bicubic_interpolation): the sixteen surrounding points are used. 

##### Step 3: Thresholding
The negative value of the rotated image correspond to the foreground. The rest is considered the background. The resulting binary image is then saved.

### Compile and run
##### Installation of DGtal and additional libraries
Concerning  the 2D part of the project, no particular install option is required. 
For the 3D part, the following libraries are needed :

 - OpenGL
 - QT4
 - QGLViewer for QT4

Additionally, it is necessary to modify the installation of DGtal. If DGtal was already installed prior to using this program, a clean uninstall might necessary. 
In the `cmake/`directory of DGtal, and for the  `FindQGLVIEWER.cmake` file :

- locate the `find_library(QGLVIEWER_LIBRARY_RELEASE ...... )` line
- add `QGLViewer-qt4` to the list of `NAMES` :
<pre><code>
find_library(QGLVIEWER_LIBRARY_RELEASE   
  NAMES qglviewer-qt4 qglviewer QGLViewer QGLViewer2 <b>QGLViewer-qt4</b>
  PATHS 
  /usr/lib
  /usr/local/lib
  /Library/Frameworks/
  ENV QGLVIEWERROOT
  ENV LD_LIBRARY_PATH
  ENV LIBRARY_PATH
  PATH_SUFFIXES QGLViewer QGLViewer/release
  )
</code></pre>
DGtal must then be compiled with QGLViewer. To do so, create the makefile with the following command:
`cmake .. -DWITH_QGLVIEWER=true`
It is also possible to turn the option on graphically using ccmake.
The installation is done the usual way:
```
make
sudo make install
```
The 3D part of the program is now runnable.

##### Compilation steps 

From the root directory of the program: 
```bash
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
The program will store the output files in `output/` file in the `rotation2D/` directory (if the program was ran from the `build/` directory. 

The `pre_processing/` folder allows to visualise the steps necessary to the preparation of the data for interpolation.
The `rotation_NN/`, `rotation_BIL/` and `rotation_BIC/` folders contain respectively the results of the rotation using nearest neighbor interpolation, bilinear interpolation and bicubic interpolation, in the .pgm and .eps formats.




