/*
 * m_options.cc
 *
 *  Created on: Jun 5, 2016
 *      Author: hazel
 */



class PySimOptions::PySimOptions: class SimOptions{


}



void PySimOptions::getSimulationMode(void) =   getLongAttr(system_options, simulation_mode, &simulation_mode );


class PySimOptions: public SimOptions {
public:
	void calculateEnergy(void);


	void getSimulationMode(void) =   getLongAttr(system_options, simulation_mode, &simulation_mode );
	void functionTwo(void);

protected:
	PyObject python_settings;
};

#endif /* SYSTEM_SIMOPTIONS_H_ */
