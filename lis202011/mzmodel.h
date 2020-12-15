#pragma once
//#include <stdio.h>
#include "consts.h"
#include <stdlib.h>
using namespace std;

enum PointStatus {
	Init = 0,
	Over_Flow = 3,
	Weir_Flow = 1,
	Orifice_Flow = 2
};

enum DischargeRuleType {
	UniformRules = 0,
	CustomRules = 1
};

class DischargeRules {
public:
	int rulesCount = 0;
	double* upLimit = NULL;
	double* lowLimit = NULL;
	PointStatus* status = NULL;
	double* factors = NULL;
};

class  CouplePoint {
public:
	//string id;//Node name
	int valid ;
	int nIndex;//Node Index
	int xIndex;//DEM col
	int yIndex;//DEM row
	int zIndex = 0;//floor
	double floodDepth;
	double overflowVol;
	double weirB;
	double cqA;
	double pondedA;
	PointStatus status;
	int factorIndex;
	int* rules;//´ýÀ©Õ¹
	double* factors;
	double oldFlow;
	double newFlow;

public:
	void updateStatus(double lastTstep, double overflow = 0.0);
	void updateDischarge(double tstep);
	void updatePointDepth(double tstep);
	void setInflow();

public:
	DischargeRules static uniformRules;
	int static loadUniformRules(string path);
	int static count;
};