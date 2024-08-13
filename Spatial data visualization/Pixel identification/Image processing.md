Useful pixels are identified by dividing the tissue scan image into 50x50 or 100x100 pixel squares that align with Patho-DBiT microfluidic barcoded pixels. The intensity within each square is calculated, and only those with signals above a set threshold are retained.

The Photoshop software will be used to generate a black and white image for pixel identification following the following steps:

1. Use Photoshop to crop the tissue scan image to maintain only region of interest that Patho-DBiT covered. The upperleft of the image should be the 1x1 pixel, and the lowerright is the 50x50. No space is allowed. See Figure 1 for example.

2. Use threashold function under Image->adjustment menu to adjust the image, so that your tissue is black and background is compeletely white.
3. Invert the color of the image. The final image is like "FFPE-2_BW.jpg" in the Example_Data folder.
4. Run the matlab script and a "postion.txt" file will be generated, which contains only the useful pixels.
   
![LM0623-Large-No barcoding-ROI](https://github.com/user-attachments/assets/10af3fa4-0220-45b2-8173-42d3167b311e)
Figure 1. Tissue scan of the region of interest

![LM0623-Large-BW](https://github.com/user-attachments/assets/00bbc08b-c09d-44bd-971d-d067eac1665b)
Figure 2. Black and white version of the tissue scan
