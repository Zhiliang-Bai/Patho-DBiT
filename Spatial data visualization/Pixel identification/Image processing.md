Useful pixels are identified by dividing the tissue scan image into 50x50 or 100x100 pixel squares that align with Patho-DBiT microfluidic barcoded pixels. The intensity within each square is calculated, and only those with signals above a set threshold are retained.

Follow these steps to generate a black and white image for pixel identification using Photoshop:

1. Crop the image:
Open the tissue scan image in Photoshop and crop it to retain only the region of interest covered by Patho-DBiT. Ensure that the upper left corner of the image corresponds to the 1x1 pixel and the lower right to the 50x50 pixel. There should be no extra space. Refer to Figure 1 for an example.

2. Adjust image threshold:
Navigate to the Image -> Adjustments -> Threshold menu to adjust the image. Modify the threshold so that the tissue area appears black and the background is completely white.

3. Invert the image colors:
Invert the colors of the image, resulting in the tissue area becoming white and the background black. The final image should resemble Figure 2.

4. Run MATLAB script:
Execute the MATLAB script 'Pixel_identification.m'. This will generate a "position.txt" file containing only the useful pixels identified from the image.
   
![LM0623-Large-No barcoding-ROI](https://github.com/user-attachments/assets/10af3fa4-0220-45b2-8173-42d3167b311e)
Figure 1. Tissue scan of the region of interest

![LM0623-Large-BW](https://github.com/user-attachments/assets/00bbc08b-c09d-44bd-971d-d067eac1665b)
Figure 2. Black and white version of the tissue scan
