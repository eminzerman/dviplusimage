# sim2convert

This shader code is the Matlab implementation of the shader presented in Easy HDR Player software, provided with SIM2 display. With the help of this code, the HDR images can be converted to a pixel-mapped LogLuv image and can be viewed on SIM2 displays using the built-in rendering mode. 

This software can be used to convert HDR images and create the SDR LogLuv coded images before subjective experiments or demonstrations.

## Use

After downloading **sim2convert.m** file, pass the HDR image to the function, e.g. **sim2convert(hdrIm)**. The output image can be saved as a bitmap of PNG image file and can be displayed on a SIM2 display using the built-in rendering mode. 

## Disclaimer

The colorspace transfer function from RGB to XYZ may be different for each SIM2 display depending on their display characteristics or the color calibration setting of the host/server computer. Therefore, you may need to change the **_mRGB2XYZ_** matrix within the code to obtain best results.

### Licence

This software is provided under GNU General Public Licence.
