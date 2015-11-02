Level 2: Photon Spectrum, Particle Evolution, Static Scenario
=============================================================

Now, we will adapt the model from level 1 to allow for a particle evolution.
That means that the particle distribution which generates the photon spectrum
is not provided by an external file, but generated from a physical model.

In this section, this physical model is as follows:

 - an astrophysical accelerator is injecting particles at a constant rate 
   distributed in a power law with constant index.

 - the ambient medium (magn. field strength, density etc.) is constant in time 

To calculate the particle evolution, we will use the ``Particles`` class.
The following example will be displayed in ``python`` style but of course
the ``C++`` syntax is analogue - see level 1 as an example.

You can find the python script ``TutorialLevel2.py`` and parameter file 
``TutorialLevel3Params.dat`` in the``docu/tutorial/`` directory. 
Run it like::

  $ python docu/tutorial/TutorialLevel2.py docu/tutorial/TutorialLevel2Params.dat

This will create a plot ``docu/pics/TutorialLevel2.png``


Let's go through the code:

 - The parameter file is read in using the ``ConfigParser`` module.

  .. sourcecode:: python

    import ConfigParser

    if __name__ == "__main__":

      # Read in parameter file
      configFile = os.path.abspath(sys.argv[1])
      configParser = ConfigParser.RawConfigParser()
      configParser.read(configFile)
      lum = float(configParser.get('Parameters','Luminosity'))
      age = float(configParser.get('Parameters','Age'))

  and so forth.

 - intialise objects:

  .. sourcecode:: python

      fr = gamerapy.Radiation()
      fp = gamerapy.Particles()
      fu = gamerapy.Utils()
      fu.DrawGamera()

  Yes, you should draw Gamera on your terminal!


 - set parameters:

  .. sourcecode:: python

      fp.SetLuminosity(lum)
      fp.SetBField(bfield)
      fp.SetEmax(emax)
      fp.SetEmin(emin)
      fp.SetSpectralIndex(spind)
      fp.SetEnergyBinsForNumericalSolver(ebins)
      fp.SetAmbientDensity(dens)
      fp.SetAge(age) 

      fr.SetDistance(dist)
      fr.SetAmbientDensity(fp.GetAmbientDensity())
      fr.SetBField(fp.GetBField())
      fr.AddThermalTargetPhotons(2.7,0.25*1.602*1.e-12)
      fr.AddThermalTargetPhotons(t,e)

  .. note::
    
    The input units for ``GAMERA`` are **always**:

    - Energy: erg
    - Time: yrs
    - Length: cm
    - Mass: g
    - Density: cm^-3
    - B-field: G


  .. hint::
    Setting the energy bins of the particle spectrum via 
    ``SetEnergyBinsForNumericalSolver`` is optional. The default number of bins
    is 100. The less bins the faster the computation, but also the less precise 
    the result! Vice versa, if the particle spectrum computation is not a time-
    critical step, the higher the better.

 - calculate a lookup table that holds the electron cooling rate due to the 
   IC process and give it to the ``Particle`` objects

  .. sourcecode:: python
    
      fr.CreateICLossLookup()
      fp.SetICLossLookup(fr.GetICLossLookup())


 - compute the particle evolution and give the result to the ``Radiation`` object:

  .. sourcecode:: python

    fp.CalculateElectronSpectrum()
    fr.SetElectrons(fp.GetParticleSpectrum())

 - with the last line, the electron distribution was computed. Doing the same 
   for protons (just for fun with somewhat different parameters) is easy:

  .. sourcecode:: python
    
    fp.SetLuminosity(1.e2*lum)
    fp.SetEmax(1.e3)
    fp.SetSpectralIndex(spind-0.1)
    fp.CalculateProtonSpectrum()
    fr.SetProtons(fp.GetParticleSpectrum())

 - now calculate the radiation spectra the provided particles emit with

  .. sourcecode:: python

    fr.CalculateDifferentialPhotonSpectrum()

 - and get the result

  .. sourcecode:: python

    TotalSED = np.array(fr.GetTotalSED())
    ICSED = np.array(fr.GetICSED())
    BremsSED = np.array(fr.GetBremsstrahlungSED())
    SynchSED = np.array(fr.GetSynchrotronSED())
    PPSED = np.array(fr.GetPPSED())
    ElectronSED = np.array(fp.GetParticleSED())

  .. hint::
    Input and output in gamerapy are ``doubles`` and 
    two-dimensional ``lists`` or `numpy <http://www.numpy.org/>`_ ``arrays``.
  
Here is the output, plotted with `matplotlib <http://matplotlib.org/>`_ :


.. image::  ../../pics/TutorialLevel2.jpg
   :align:  center
   :scale: 100 %


