/* -*- mode:c++ -*- */
#ifndef BASE_MANAGER_H
#define BASE_MANAGER_H

//#include "Utility.h"
#include "Parameters.h"
#include "Data.h"
//#include "Kernel.h"

using namespace std;

/**
 * @class BaseManager
 *
 * @brief Base class for all action managers.
 * Provides a container for the Data and the Parameters classes.
 * Provides output facilities.
 */
class BaseManager {
protected:
	Parameters* mpParameters;
	Data* mpData;
public:
	BaseManager(Parameters* apParameters, Data* apData);
	void Init(Parameters* apParameters, Data* apData);
};

#endif /* BASE_MANAGER_H */
