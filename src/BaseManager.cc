#include "BaseManager.h"

//BaseManager::BaseManager(Parameters* apParameters, Data* apData) {
//	Init(apParameters, apData);
//}
//
//void BaseManager::Init(Parameters* apParameters, Data* apData) {
//	mpParameters = apParameters;
//	mpData = apData;
//}

BaseManager::BaseManager(Parameters* apParameters) {
	Init(apParameters);
}

void BaseManager::Init(Parameters* apParameters) {
	mpParameters = apParameters;
}
