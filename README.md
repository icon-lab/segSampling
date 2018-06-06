# segSampling

This is a collection of Matlab functions and demos to demonstrate the techniques described in the following paper:

Senel LK, Kilic T, Gungor A, Kopanoglu E, Guven HE, Saritas EU, Koc A, Çukur T. Statistically segregated k-space sampling for accelerating multiple-acquisition MRI. arXiv:1710.00532, 2017.

You are encouraged to modify/distribute this code. However, please acknowledge this code and cite the paper appropriately.
This code utilizes and requires SPIRiT library.


## Demos
```matlab
	demo_segsampling_in_vivo_T2_PE.m : Demonstrates segregated sampling on an in vivo T2-weighted dataset.
        demo_segsampling_in_vivo_bSSFP_PE.m : Demonstrates segregated sampling on in vivo phase-cycled bSSFP datasets.
	demo_segsampling_phantom_T1_ZF.m : Demonstrates segregated sampling on a phantom T1-weighted dataset.
	demo_segsampling_phantom_T2_PE.m : Demonstrates segregated sampling on phantom phase-cycled bSSFP datasets.

   note: Please run setPath.m prior to the demo scripts.
```

For any questions, comments and contributions, please contact
Tolga Cukur (cukur[at]ee.bilkent.edu.tr)

(c) ICON Lab 2017
