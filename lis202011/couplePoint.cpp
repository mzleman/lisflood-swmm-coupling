#include "global.h"
#include "couplefuncs.h"



using namespace std;

void noPondingOverflow()
{
	int i;
	double overflow;
	for (i = 0; i < CouplePoint::count; i++) {
		overflow = getOverflow(QPoints[i].nIndex);
		QPoints[i].updateStatus(overflow);
	}
}

//void pondingOverflow() //！暂未使用
//{
//	int i;
//	for (i = 0; i < CouplePoint::count; i++)
//	{
//		QPoints[i].updateStatus();//不传入参数
//	}
//
//}


void CouplePoint::updateStatus(double lastTstep, double overflow)
{
	if (this->valid && isEmptyDEMGrid(this->yIndex, this->xIndex)) {
		cout << "点 [" << this->yIndex << "," << this->xIndex << "] 不合法" << endl;
		this->valid = FALSE; // 耦合点无高程值, 标记为无效耦合点
		return;
	}
	if (!Parptr->allowPonding) //服务于noponding overflow，不允许雨水井积水时，使用swmm计算的overflow流量，更新溢流总体积。
	{
		double floodDepth;

		if (overflow > 0.0)
		{
			this->overflowVol = overflow * lastTstep; //计算溢流量
		}
		else
		{
			this->overflowVol = 0.0;
		}

		this->floodDepth = Arrptr->H[this->xIndex + this->yIndex * Parptr->xsz];

	}
	else   //服务于ponding overflow，允许雨水井积水时，更新水头差
	{
		this->overflowVol = 0.0;
		int overFlowFlag;
		double head = getNodeHead(this->nIndex, &overFlowFlag);
		if (overFlowFlag)
		{
			this->floodDepth = Arrptr->H[this->xIndex + this->yIndex * Parptr->xsz] +
				Arrptr->DEM[this->xIndex + this->yIndex * Parptr->xsz] - head;                    // 地表水位=水深+地表高程， 水头差=地表水位-管点积水水位
		}
		else
		{
			this->floodDepth = Arrptr->H[this->xIndex + this->yIndex * Parptr->xsz];         //管点不溢流，则水头差为地表积水水深，而非减去管点水位。
		}
	}
}


int getRulesIndex(double depth)
{
	int i, index = -1;
	for (i = 0; i < CouplePoint::uniformRules.rulesCount; i++) {
		if (depth >= CouplePoint::uniformRules.lowLimit[i] && depth < CouplePoint::uniformRules.upLimit[i]) {
			index = i; break;
		}
	}
	return index;
}


void CouplePoint::updateDischarge(double tstep) //更新耦合点排涝流量或者溢流量
{
	if (!this->valid) return;
	double floodDepth, maxVol;
	int ruleIndex;
	int default_flag;
	double g = 9.80;
	//if ( isEmptyDEMGrid(this->yIndex, this->xIndex) ) //耦合点无高程值
	//{
	//	return;
	//}
	this->oldFlow = this->newFlow;
	if (!Parptr->allowPonding && (this->overflowVol > EPS)) //overflow
	{
		this->newFlow = 0.0;
	}
	else if (this->floodDepth < 0.0)// 注意,allowPonding==OFF时,floodDepth一定大于等于0。所以allowPonding==OFF&& overflowVol==0时,只会走discharge
	{
		/*default_flag = 1;*/
		/*this->newFlow = 0.0;*/
		floodDepth = -1.0*this->floodDepth;
		this->overflowVol = floodDepth * this->pondedA;
		this->newFlow = this->overflowVol / tstep;
	}
	else //discharge
	{
		default_flag = TRUE;
		floodDepth = this->floodDepth;
		maxVol = floodDepth * Parptr->dx * Parptr->dy;
		if (Parptr->dischargeRules == UniformRules) //uniform rules 所有节点使用统一的流态排水规则
		{
			ruleIndex = getRulesIndex(floodDepth);
			if (ruleIndex == -1); //not found rule, do nothing here,use default,
			else
			{
				switch (CouplePoint::uniformRules.status[ruleIndex])
				{
				case Weir_Flow:
					this->newFlow = (-1.0)*CouplePoint::uniformRules.factors[ruleIndex] * this->weirB*pow(floodDepth, 1.5)*pow(2 * g, 0.5);
					break;
				case Orifice_Flow:
					this->newFlow = (-1.0)*CouplePoint::uniformRules.factors[ruleIndex] * this->cqA*pow(2 * g * floodDepth, 0.5);
					break;
				default:
					this->newFlow = 0.0;
					break;
				}
				default_flag = FALSE;//not use default
			}
		}
		else if (Parptr->dischargeRules == CustomRules)//custom rules
		{
			/*待扩展*/
			/*待扩展*/
			default_flag = FALSE;//not use default
		}
		if (default_flag) //default rules
		{
			/*floodDepth = domain.waterDepths.cellArray[this->yIndex][this->xIndex].values[this->zIndex];*/
			floodDepth = this->floodDepth;
			maxVol = floodDepth * Parptr->dx * Parptr->dy;
			if (floodDepth > 0 && floodDepth <= DEFAULT_THRESHOLD)
			{
				this->newFlow = (-1.0)*DEFAULT_WEIR_COEFF * this->weirB*pow(floodDepth, 1.5)*pow(2 * g, 0.5);
			}
			if (floodDepth > DEFAULT_THRESHOLD)
			{
				this->newFlow = (-1.0)*DEFAULT_CQ * this->cqA*pow(2 * g * floodDepth, 0.5);
			}
		}
		if (maxVol < this->newFlow * tstep) //防止排水量大于栅格内水量
		{
			this->newFlow = (-1.0)*maxVol / tstep;
		}
	}
}

void CouplePoint::updatePointDepth(double tstep) //排水,更新地面水深
{
	if (!this->valid) return;
	double A = Parptr->dx * Parptr->dy;
	if (this->overflowVol > 0.0) //overflow
	{
		Arrptr->H[this->xIndex + this->yIndex * Parptr->xsz] += this->overflowVol / A;
		//this->overflowVol = 0.0;                    //重置为0.0
	}
	else
	{
		Arrptr->H[this->xIndex + this->yIndex * Parptr->xsz] += this->newFlow * tstep / A;
		/*wd[this->zIndex] += ( (this->newFlow + this->oldFlow)/2.0*modelControl.timeStep/a );*/
	}
}


void CouplePoint::setInflow() // 与SWMM交互，为雨水点注入排涝流量
{
	if (!this->valid) return;
	setLatFlow(this->nIndex, -1.0*this->newFlow);
}


int CouplePoint::loadUniformRules(string path) // 加载全局排涝计算规则
{
	ifstream rulesFile(path);
	string temp;
	int count = 0;
	if (!rulesFile) {
		cout << "uniformRules File open fail" << endl;
		return -1;
	}
	rulesFile >> temp >> count;
	CouplePoint::uniformRules.rulesCount = count;
	/*do {
		getline(rulesFile, temp);
	} while (*temp.c_str() != ';');*/

	double* upls =
		CouplePoint::uniformRules.upLimit = new double[count];

	double* lowls =
		CouplePoint::uniformRules.lowLimit = new double[count];

	PointStatus* pss =
		CouplePoint::uniformRules.status = new PointStatus[count];

	double* fcts =
		CouplePoint::uniformRules.factors = new double[count];

	//读取rules
	cout << "\nUniform Rules:" << endl;
	for (int i = 0; i < count; i++)
	{

		rulesFile >> lowls[i] >> upls[i] >> temp >> fcts[i]; //读取一行记录

		if (!temp.compare("WEIR")) pss[i] = Weir_Flow;
		else if (!temp.compare("ORIFICE")) pss[i] = Orifice_Flow;
		else pss[i] = Init;
		cout << "\tlow:" << lowls[i] << "\tup:" << upls[i] << "\tstatus:" << pss[i] << "\tfactor:" << fcts[i] << endl;
	}
	rulesFile.close();
	return 0;
}