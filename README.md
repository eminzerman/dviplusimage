# dviplusimage

This project consists of the Matlab code of the HDR image rendering algorithm. This code is developed to render HDR image using SIM2 displays' DVI PLUS mode, which enables users to present HDR images with custom rendering. 

This software can be used as a tool for subjective experiments and psychometric studies, as it reproduces HDR images accurately considering the parameters of the display. 

## Use

You can simply clone or download the project and run **dviplusimage** within Matlab. Details are presented in each m-file. There are three options:
* **dviplusim.m** file is a single-file Matlab function which can be used to render HDR images. It only requires two files: *led_psf.mat* and *lcdGammaCorr.mat*.
* **dviplusvid.m** file is a single-file Matlab function which can be used to render HDR images.
* **dviplusimage** folder includes a multpile-file version of the same code. **dviplusimage\dviplusimage.m** can be used the same way as the fist option for natural images. Alternatively, **dviplusimage\dviplusimageTest.m** should be used for *test* images with sharp discontinuities (i.e. luminance masking, contrast studies).

## Disclaimer

Specifically, this particular version is developed considering SIM2 HDR47 E S4K model. The software was validated using only one SIM2 display. Therefore, special care may need to be taken for color calibration before using the software. 

Although it is not tested, the rendering algorithm can also be used for other HDR displays, considering that some key parameters (e.g. maximum luminance, point spread function, hardware power limitations, gamma functions for each color channel, etc.) are provided.

### Licence

This software is provided under GNU General Public Licence.

If you use this software for research purposes, please kindly cite our papers below:
* E. Zerman, G. Valenzise, F. De Simone, F. Banterle, and F. Dufaux. **_"Effects of Display Rendering on HDR Image Quality Assessment."_** SPIE Optical Engineering+ Applications, Applications of Digital Image Processing XXXVIII, San Diego, CA, USA, August 2015.
* E. Zerman, G. Valenzise, and F. Dufaux. **_"A Dual Modulation Algorithm for Accurate Reproduction of High Dynamic Range Video."_** IEEE 12th Image, Video, and Multidimensional Signal Processing Workshop (IVMSP), Bordeaux, France, July 2016.