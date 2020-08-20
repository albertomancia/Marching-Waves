# Marching-Waves

Copyright © 2020 Alberto Mancia. All rights reserved.

Art program that utilizes James Sethian's fast marching method for simulating wave propagation using the Eikonal Equation.

![Sample output](https://raw.githubusercontent.com/albertomancia/Marching-Waves/master/example.png)

# How to use

Requires Processing 3.5.4: http://processing.org

Load an image by putting it in the same folder as the program, and changing the string called "filename" to match. The image you load into the program must satisfy the following:
  1. The background, or any part of the image that you don't want to draw on, must be #0000FF (0 red, 0 green, 255 blue)
  2. Any pixels that you want to use as the starting set must be #FF0000 (255 red, 0 green, 0 blue). This doesn't need to be included, as you can draw your own origin set with the cursor.
  3. The rest of the image should be grayscale to ensure that the program doesn't misread it.
  
When the program starts, it is in drawing mode. You can add to the origin set by clicking and dragging with your mouse, or by pressing one of the following keys:

* 'r' to place random dots across the canvas
* 'l' to draw random lines across the canvas
* 'g' to draw a square grid across the canvas (change the "gridsize" variable to adjust)

Once you've drawn the origin set, press ENTER to begin calculation. After the calculation finishes, by default the program will show the animated lines moving out from the origin set. At any time, you can press ENTER again to reset the program and start over.

In Line View, press SPACE to play/pause the animation and use the left and right arrow keys to move between frames.

To view the vector and scalar fields generated by the algorithm, press 'f' to enter Field View and use 'q' to cycle between:

* The solution to the Eikonal Equation, as a 2D scalar field. This represents the heightmap that is being sliced by the contour lines, but more fundamentally it represents the path-distance from the origin set for each pixel, using the reference image as a cost field. It represents the two ways this program's output can be interpreted: as a valley whose walls are made steeper or shallower by the image beneath, or as a series of waves whose speed is changed by the image as they pass over it.
* The gradient field of the solution, or the uphill direction of the valley. This is represented by mapping the direction to a hue along the color wheel.
* The curl of the gradient field, or how much the gradient vectors are changing direction around a given point.

# Exporting SVGs (beta)

In Line View, pause the animation and select a frame using the arrow keys, then press 's' to save an SVG of the frame. This works using the fields described above, dithering along the lines using the curl field to distribute anchor points , and then using the tangent vector (90° rotation of the gradient vector) at each point to create the appropriate Bezier vertex.

Currently working on a new method that uses marching squares to draw contour lines along the solution field, and group the line segments using a modified bubble sort algorithm.
