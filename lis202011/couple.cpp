/*
	load Swmm Dll
	load coupling functions
*/

#include "couplefuncs.h"

using namespace std;

//SWMM functions
HINSTANCE SwmmDll = NULL;
GetCouplePointsN getCouplePointsN = NULL;
GetCouplePoints getCouplePoints = NULL;
SWMM_OPEN swmmOpen = NULL;
SWMM_START swmmStart = NULL;
GetSWMMTotalTime getSWMMTotalTime = NULL;
SWMM_STEP swmmStep = NULL;
SWMM_END swmmEnd = NULL;
SWMM_CLOSE swmmClose = NULL;
GetOverflow getOverflow = NULL;
SetLatFlow setLatFlow = NULL;
GetSWMMTstep getSWMMTstep = NULL;
SetAllowPonding setAllowPonding = NULL;
GetNodeHead getNodeHead = NULL;
ReportNodeFlood reportNodeFlood = NULL;

//declare relevant variables
DischargeRules  CouplePoint::uniformRules;
CouplePoint* QPoints = NULL;
int CouplePoint::count = 0;

int loadSWMM() {
	SwmmDll = LoadLibrary(TEXT("SwmmSourceCode.dll"));
	swmmOpen = (SWMM_OPEN)GetProcAddress(SwmmDll, "swmm_open");
	swmmStart = (SWMM_START)GetProcAddress(SwmmDll, "swmm_start");
	swmmStep = (SWMM_STEP)GetProcAddress(SwmmDll, "swmm_step");
	swmmEnd = (SWMM_END)GetProcAddress(SwmmDll, "swmm_end");
	swmmClose = (SWMM_CLOSE)GetProcAddress(SwmmDll, "swmm_close");
	getCouplePointsN = (GetCouplePointsN)GetProcAddress(SwmmDll, "swmm_getCouplePointsN");
	getCouplePoints = (GetCouplePoints)GetProcAddress(SwmmDll, "swmm_getCouplePoints");
	getSWMMTotalTime = (GetSWMMTotalTime)GetProcAddress(SwmmDll, "swmm_getSWMMSimTime");
	getOverflow = (GetOverflow)GetProcAddress(SwmmDll, "swmm_getOverflow");
	setLatFlow = (SetLatFlow)GetProcAddress(SwmmDll, "swmm_setLatFlow");
	getSWMMTstep = (GetSWMMTstep)GetProcAddress(SwmmDll, "routing_getRoutingStep");
	setAllowPonding = (SetAllowPonding)GetProcAddress(SwmmDll, "swmm_setOption_allowPonding");
	getNodeHead = (GetNodeHead)GetProcAddress(SwmmDll, "swmm_getNodeHead");
	reportNodeFlood = (ReportNodeFlood)GetProcAddress(SwmmDll, "swmm_nodeFlood");

	if (!SwmmDll) {
		cout << "load swmm fail" << endl;
		return -1;
	}
	if (!swmmOpen) {
		cout << "load swmm_open fail" << endl;
		return -1;
	}
	if (!swmmStart) {
		cout << "load swmm_start fail\n" << endl;
		return -1;
	}
	if (!swmmStep) {
		cout << "load swmm_step fail\n" << endl;
		return -1;
	}
	if (!swmmClose) {
		cout << "load swmm_close fail\n" << endl;
		return -1;
	}
	if (!swmmEnd) {
		cout << "load swmm_end fail\n" << endl;
		return -1;
	}
	if (!getCouplePointsN) {
		cout << "load swmm_getCouplePointsN fail" << endl;
		return -1;
	}
	if (!getCouplePoints) {
		cout << "load swmm_getCouplePoints fail" << endl;
		return -1;
	}
	if (!getSWMMTotalTime) {
		cout << "load swmm_getSWMMSimTime fail" << endl;
		return -1;
	}
	if (!getOverflow) {
		cout << "load swmm_getOverflow fail" << endl;
		return -1;
	}
	if (!setLatFlow) {
		cout << "load swmm_setLatFlow fail" << endl;
		return -1;
	}
	if (!getSWMMTstep) {
		cout << "load routing_getRoutingStep fail" << endl;
		return -1;
	}
	if (!setAllowPonding) {
		cout << "load swmm_setOption_allowPonding fail" << endl;
		return -1;
	}
	if (!getNodeHead) {
		cout << "load swmm_getNodeHead fail" << endl;
		return -1;
	}
	if (!reportNodeFlood) {
		cout << "load reportNodeFlood fail" << endl;
		return -1;
	}

	cout << "loadSWMM success" << endl;
	return 0;
}

int initPoints() {
	int n = 0;
	/*NCouplePoint = 0;*/
	n = getCouplePointsN(Parptr->blx, Parptr->bly, Parptr->ysz, Parptr->xsz, Parptr->dx);
	int* indexs = new int[n];
	int* rows = new int[n];
	int* cols = new int[n];
	double* cqAs = new double[n];
	double* weirBs = new double[n];
	double* pondedAs = new double[n];
	/*NCouplePoint=getCouplePoints(domain.blx, domain.bly, domain.rows, domain.cols, domain.dx, indexs, rows, cols, cqAs, weirBs);*/
	CouplePoint::count = getCouplePoints(Parptr->blx, Parptr->bly, Parptr->ysz, Parptr->xsz, Parptr->dx, indexs, rows, cols, cqAs, weirBs, pondedAs);
	if (CouplePoint::count != n) {
		cout << "couple points count wrong" << endl;
		return -1;
	}
	QPoints = (CouplePoint*)calloc(n, sizeof(CouplePoint));
	if (!QPoints) {
		cout << "calloc QPoints fail" << endl;
		return -1;
	}
	for (int i = 0; i < n; i++) {
		QPoints[i].valid = TRUE;
		QPoints[i].nIndex = indexs[i];
		QPoints[i].yIndex = rows[i];
		QPoints[i].xIndex = cols[i];
		QPoints[i].zIndex = 0;
		QPoints[i].cqA = cqAs[i] > 0.0 ? cqAs[i] : DEFAULT_ORIFICE_AREA;
		QPoints[i].weirB = weirBs[i] > 0.0 ? weirBs[i] : DEFAULT_WEIR_B;
		QPoints[i].floodDepth = 0.0;
		QPoints[i].oldFlow = 0.0;
		QPoints[i].newFlow = 0.0;
		QPoints[i].overflowVol = 0.0;
		QPoints[i].status = Init;
		QPoints[i].rules = NULL;
		QPoints[i].factorIndex = -1;
		QPoints[i].factors = NULL;
		QPoints[i].pondedA = pondedAs[i];
		//cout << "节点索引: " << indexs[i] << "\t行号: " << rows[i] << "\t列号: " << cols[i] << "\tcqA: " << QPoints[i].cqA << "\tweirB: " << QPoints[i].weirB << "\tpondedA: " << QPoints[i].pondedA <<endl;
	}
	delete[] indexs; indexs = NULL;
	delete[] rows; rows = NULL;
	delete[] cols; cols = NULL;
	delete[] cqAs; cqAs = NULL;
	delete[] weirBs; weirBs = NULL;
	cout << "\ninit Points success" << endl;
	return 0;
}

int initCouple() {
	string inpPath = Fnameptr->inpFileName;
	string rptPath, outPath;
	rptPath = inpPath;
	rptPath.replace(inpPath.rfind(".inp"), 4, ".rpt");
	outPath = inpPath;
	outPath.replace(inpPath.rfind(".inp"), 4, ".out");
	if ( loadSWMM() ) {
		return -1;
	}
	if ( swmmOpen(inpPath.c_str(), rptPath.c_str(), outPath.c_str()) ) {
		cout << "inp文件打开失败" << endl;
		return -1;
	}
	if ( swmmStart(1) ) {
		cout << "swmm启动失败" << endl;
		return -1;
	}
	if ( initPoints() ) {
		cout << "初始化耦合点失败" << endl;
		return -1;
	}
	if (Parptr->dischargeRules == UniformRules) {
		CouplePoint::loadUniformRules(Fnameptr->uniformRulesFileName);
	}
	setAllowPonding(Parptr->allowPonding);
	return 0;
}

void updateTimeStep(int swmm) {
	if (Statesptr->acceleration || Statesptr->Roe) {
		if (Statesptr->acceleration == ON) CalcT(Parptr, Solverptr, Arrptr);//计算时间步长，修改Soverlptr->Tstep
		if (Statesptr->Roe == ON) CalcTRoe(Parptr, Solverptr, Arrptr);//计算时间步长
		if (swmm) {
			double swmmTstep;
			swmmTstep = getSWMMTstep(DW, Solverptr->InitTstep);
			Solverptr->Tstep = getmin(swmmTstep, Solverptr->Tstep);
		}
	}
}

void updateCouplePoints(double tstep, double lastTstep) {
	double overflow;
	#pragma omp parallel for private(overflow)
	for (int i = 0; i < CouplePoint::count; i++)
	{
		overflow = getOverflow(QPoints[i].nIndex);
		QPoints[i].updateStatus(lastTstep, overflow); // 上一次迭代的计算步长
		QPoints[i].updateDischarge(tstep); // 当前的计算步长
		QPoints[i].updatePointDepth(tstep);
		QPoints[i].setInflow();
	}
}