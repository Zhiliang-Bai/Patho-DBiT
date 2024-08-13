Useful pixels are identified by dividing the tissue scan image into 50x50 or 100x100 pixel squares that align with Patho-DBiT microfluidic barcoded pixels. The intensity within each square is calculated, and only those with signals above a set threshold are retained.

There two steps: To run the Matlab script "Pixel_identification.m"

1. Use Photoshop or other photo editing software to crop the microscope image into exactly the size of the DBiT-seq covering area. For example, the upperleft of the image should be the 1x1 pixel of DBiT-seq, and the lowerright is the 50x50. No space is allowed. See "FFPE-2.jpg" for example.

2. Use threashold function under Image->adjustment menu to adjust the image, so that your tissue is black and background is compeletely white.
3. Invert the color of the image. The final image is like "FFPE-2_BW.jpg" in the Example_Data folder.
4. Run the matlab script and a "postion.txt" file will be generated, which contains only the useful pixels.
   
![LM0623-Large-No barcoding-ROI](https://github.com/user-attachments/assets/10af3fa4-0220-45b2-8173-42d3167b311e)
Figure 1. Tissue scan of the region of interest

![LM0623-Large-BW](https://github.com/user-attachments/assets/00bbc08b-c09d-44bd-971d-d067eac1665b)
Figure 2. Black and white version of the tissue scan
