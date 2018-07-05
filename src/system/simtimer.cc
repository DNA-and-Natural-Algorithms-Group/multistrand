#include "simoptions.h"
#include "simtimer.h"


SimTimer::SimTimer(SimOptions& myOptions) {

	maxsimtime = myOptions.getMaxSimTime();
	stopcount = myOptions.getStopCount();
	stopoptions = myOptions.getStopOptions();

	// saving the pointer to enable access to cotranscriptional timing values
	simOptions = &myOptions;

}

// advances the simulation time according to the set rate
void SimTimer::advanceTime(void) {

	rchoice = rate * drand48();
	stime += (log(1. / (1.0 - drand48())) / rate);

}

// returns TRUE if this transition needs to be executed.
bool SimTimer::wouldBeHit(const double rate) {
	return rchoice < rate;
}

// returns which Nth collision needs to be used
int SimTimer::checkHitBi(const double collisionRate) {

	return (int) floor(rchoice / collisionRate);

}

// returns TRUE if this transition needs to be executed
// and subtracts the rate
bool SimTimer::checkHit(const double rate) {

	rchoice = rchoice - rate;
	return (rchoice < 0);

}

// returns TRUE if a new nucleotide is to be added to the chain
bool SimTimer::checkForNewNucleotide(void) {

	if (simOptions->cotranscriptional && stime > (nuclAdded + simOptions->initialActiveNT) * simOptions->cotranscriptional_rate) {

		nuclAdded++;
		return true;
	}

	return false;
}

std::ostream& operator<<(std::ostream& ss, SimTimer& timer) {

	ss << "rchoice ";
	ss << timer.rchoice << "  rate  ";
	ss << timer.rate << "  simTime  ";
	ss << timer.stime << "\n";

	return ss;
}
