#/usr/local/bin/python
import numpy as np
import gappa as gp
import unittest



class TestRadiationClass(unittest.TestCase):

    def grey_body(self,e,t,e_dens):
        return 15*e_dens*(gp.pi*gp.kb*t)**-4 *e**2 / (np.exp(e /(gp.kb*t))-1)

    def write_grey_body(self,e,t,e_dens):
        g = self.grey_body(e,t,e_dens)
        f = open("rad_temp.dat",'w')
        for ee,gg in list(zip(e,g)):
            f.write(str(ee/gp.eV_to_erg)+" "+str(gg*gp.eV_to_erg)+"\n")
        f.close()

    def test_AddThermalTargetPhotons(self):
        fr = gp.Radiation()
        t = 2.7
        e_dens = 0.25*gp.eV_to_erg
        fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg)
        g_1 = np.array(fr.GetTargetPhotons())
        g_2 = self.grey_body(g_1[:,0],t,e_dens)

        dg = np.sqrt(np.sum((g_2-g_1[:,1])**2)) / np.sum(g_2)
        self.assertTrue(dg < 1e-4)

    def test_AddArbitraryTargetPhotons(self):
        fr = gp.Radiation()
        t = 2.7
        e_dens = 0.25*gp.eV_to_erg
        fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg)
        g_a = np.array(fr.GetTargetPhotons(0))

        g_2 = self.grey_body(g_a[:,0],t,e_dens)

        fr.AddArbitraryTargetPhotons(list(zip(g_a[:,0],g_2)))
        g_b = np.array(fr.GetTargetPhotons(1))
        
        g_1 = np.interp(g_a[:,0],g_b[:,0],g_b[:,1])

        dg = np.sqrt(np.sum((g_2-g_1)**2)) / np.sum(g_2)
        self.assertTrue(dg < 1e-4)


    def test_ImportTargetPhotonsFromFile(self):
        fr = gp.Radiation()
        t = 2.7
        e_dens = 0.25*gp.eV_to_erg
        fr.AddThermalTargetPhotons(2.7,0.25*gp.eV_to_erg,1000)
        g_a = np.array(fr.GetTargetPhotons(0))
        self.write_grey_body(g_a[:,0],t,e_dens)
        fr.ImportTargetPhotonsFromFile("rad_temp.dat")
        g_2 = self.grey_body(g_a[:,0],t,e_dens)
#        print g_2
        g_b = np.array(fr.GetTargetPhotons(1))
        g_1 = np.interp(g_a[:,0],g_b[:,0],g_b[:,1])
        dg = np.sqrt(np.sum((g_2-g_1)**2)) / np.sum(g_2)
        self.assertTrue(dg < 1e-4)



suite = unittest.TestLoader().loadTestsFromTestCase(TestRadiationClass)
unittest.TextTestRunner(verbosity=3).run(suite)
#if __name__ == '__main__':
#    unittest.main()
