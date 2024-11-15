11/2024

Welcome to the Wisconsin version of the "Fischer group+" receiver function and CCP stack code.
Contributions to code from: Eva Golos, Junlin Hua, Emily Hopper, Vedran Lekic, David Abt, and others…

All steps are run from the script master.m in the directory NewFunctions. I recommend running it a section at a time. The structure of this code is such that that Auto_Prep.m and Generate_RFs.m contain many functions which handle the majority of the workflow, including data downloading*, arrival picking (necessary to select analysis windows in regions where the arrival differs from predictions using 1D models via TauP), free-surface transformation, and receiver function generation and stacking. Different flags are used to call different functionality from master.m. CCP stacking is performed most currently through function ccp_stack_slopeangle.

*Recently, there are major issues with the first few sections as the IRIS event and waveform request tools have been changed (i.e. Breq_Fast is no longer supported by Earthscope). I have a few hasty fixes to get SAC files in the format we need from existing tools, but for now I recommend sticking with the Alaska_test data which is already in SAC format with the necessary metadata. I will push an update to the code in the near future which works otherwise.

New features: function Generate_RFs_RTZ.m keeps the receiver functions in R, T system - no free surface transformation.



##############
Directions:

In order to fit within limits on file sizes and number of files for github uploads, some directories are stored as zip files that need to be unzipped. These are:
Functions/m_map.zip
NewFunctions/PROPMAT.zip

The second thing you will need to do is download the seizmo package and place it in the top directory:
https://epsc.wustl.edu/~ggeuler/codes/m/seizmo/

Similarly, download the TauP package and place it in the Functions directory: 
https://www.seis.sc.edu/taup/

Now you're ready!


The steps in master.m (some of this is repeated in the comments in the script)

1. Set up paths
2. Request data (doesn't work currently)
3. Unpack data, sort into SAC files. Also remove repeat data and stations with no data
4. Create structure SAC_Filenames.mat, which is used to locate data in following steps
5. Optional: make sure events are up to date
6. Run the array picker - the most time-consuming step. Calculates predicted arrival time using TauP, but then more carefully compares waveforms from nearby stations to get a more accurate arrival pick. Also calculates signal to noise ratio
7. Reformatting of array pick data (very quick)
8. Prep data structures for receiver functions. This step also is where the free-surface transform is applied to do P-SV conversion
9. Set variables for generating receiver functions
10. Generating receiver functions. If you don't have the PROPMAT executables compiled, this will fail at the stage of "calculating reference synthetics," which are needed to do full quality control. But you will still be able to look at the raw receiver functions and do your own ad hoc QC!
11. Using the receiver functions, stack them by station and display plots. Option to sort by backazimuth, but this won't work unless you have many events.
12. Common Conversion Point stacking of RFs**. This will work even if you couldn't get PROPMAT working. Very time consuming.




##############




Some notes and tips for running:
- The first steps are related to requesting, unpacking, and prepping waveform data. Different functions within Auto_Prep.m are called. CAUTION the IRIS station search tool used here is deprecated. The code for this part needs to be changed. An example of the format the code currently works with is included as Data/Station_Lists/rawstas_Alaska_test.txt. The sections that won't work are commented out currently.

- NEIC catalogue (events) was updated December 2023. To include more recent events, set NEIC flag to ‘yes’ and run the appropriate section. Code will update and generate a new version, and save old catalogue in Data/Event_Lists/Accessory.

- An example dataset (data from the MOOS network) is included in Data/Projects/Alaska_test. Follow that to see the directory structure needed.

- The Array Picker step easily takes the most time -- maybe several days depending on how much data you have. Make sure you're running it on a machine that can be left running for a while.

- Junlin and I found we had to change the phase window in Set_Anaylsis_Windows within Auto_Prep, as well as lines 2695-2696 in Generate_RFs. For a new dataset/period band you should investigate if you need to adjust also. Function plotWindow.m can be used to visualize if the window captures the entirety of the envelope peak (run from within AutoPrep, set if clause to 1, line 3686. I’d suggest doing this in debug mode).

- There are a lot of quality-control steps that happen within Generate_RFs that are specific to mantle structure. I haven't explored these for the Ps case heavily, and they may need to be adjusted (especially for basin settings). I'd suggest starting by not using this functionality, i.e. set flag iffiltdat in Generate_RFs to 0.

- Something you might need to change is the parameters for the deconvolution in the call to function IDRF, lines ~3057 in Generate_RFs.m. For instance, changing number of iterations or Gaussian window length may impact results.

- Other software that is included but you may need to compile/check:
	- PROPMAT to calculate synthetics (lots of trouble compiling on new Mac processors; don't need until further in workflow though). This version included works for Macs pre-M1 processor (2019). I have also built it on the Rocky Linux system, please contact me for instructions on how to re-compile.
	- TauP - to calculate ray paths.
	- seizmo - though I didn't find I needed to compile anything to use.

- If you find that you're missing something (if Matlab complains that it can't find a function or variable), please let me know! Similarly, if Matlab asks you to install a toolbox, let me know so I can keep track of what's needed and improve this document. Thanks!




**A note from Junlin about CCP stacking step:
"The CCP stacking is performed with CCP_stack_slopeangle.m which use a recently developed weighting function defined from kernel shapes, and is preferred. However, for crustal structure, the weighting function does not expand horizontally much for CCP stacking with sparse array. If the station is not dense enough, the code CCP_stack.m can be used whose weighting is designed more empirically and can be adjusted. The quality control criteria used for mantle structure is also included in the code."




Good luck!

Eva Golos
University of Wisconsin - Madison, Department of Geoscience
golos@wisc.edu
