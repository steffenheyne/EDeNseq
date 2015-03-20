#include "BaseManager.h"

BaseManager::BaseManager(Parameters* apParameters, Data* apData) {
	Init(apParameters, apData);
}

void BaseManager::Init(Parameters* apParameters, Data* apData) {
	mpParameters = apParameters;
	mpData = apData;
}
