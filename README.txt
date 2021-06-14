Hi, below is information on the files in this folder, the inputs to these
functions, and how to run BAMBI!

BAMBI.m - The main function. Inputs are:
  DEMfile - an ascii file of the bathymetric grid
  flow - the flow direction
  flowvar - the variation in flow that defines a leeside. set at 40 in default.
  export - 1 is yes and 0 is no

BedformThreshold.m - a function within BAMBI that defines the scale distinction
                    between large and small scale dunes.

findWL.m - the function that measures wavelength in BAMBI.

keep.m - a function written by Xiaoning (David) Yang to keep variables defined.

slopeCALC.m - a function to calculate the slope and aspect maps from a
              bathymetric grid.

runfromexcel.m - a script to run multiple bathymetric files in batches.

parameterfile.xlsx - an example of the parameterfile input file for
                    runfromexcel.m


To run a single bathymetric map through BAMBI, you may download this folder and
save in your Documents/MATLAB folder. Then, type this line into matlab command:

[smalldunes,largedunes] = BAMBI(DEMfile,flow,flowvar,export)

If you would like to run a batch of bathymetric maps through BAMBI, fill out the
parameter file with the information required where each row represents a new file.
Then run the 'runfromexcel.m' script, which will initate BAMBI.m

If you have any questions feel free to reach out to me at jcisneros1024@gmail.com.
I am also happy to work together to improve or fix any errors. I am also very
enthusiastic to add more functionality! So feel free to commit changes, updates,
or send your ideas my way! Thanks so much!

Julia Cisneros, PhD (she/her)
NSF Postdoctoral Fellow
juliacisneros.com
