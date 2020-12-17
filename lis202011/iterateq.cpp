﻿/*
*****************************************************************************
ITERATEQ and UPDATEH
---------------------


*****************************************************************************
*/

#include "lisflood.h"
#include "global.h"
#include "update.h"

//-----------------------------------------------------------------------------
// ITERATE THROUGH TIME STEPS
void IterateQ()
{

  // Fnames *Fnameptr, Files *Fptr, States *Statesptr, Pars *Parptr, Solver *Solverptr, BoundCs *BCptr, Stage *Locptr, ChannelSegmentType *ChannelSegments, Arrays *Arrptr, vector<int> *RiversIndexVecPtr, int *RiversIndexPtr, vector<ChannelSegmentType> *ChannelSegmentsVecPtr, int *verbose
  // int i,j,ptr;
  // double FloodArea;
  // double loss; //temp variable to keep track of losses since last mass interval
  // double Comp_time, Model_Comp_Ratio, Model_time_left, Est_Time_Tot, Est_Time_Fin;


  // main iteration loop
  while (
    // we continue if
    // if simulation time is not reached
    Solverptr->t < Solverptr->Sim_Time
    // if wallclock time is not reached
    && !(Statesptr->killsim==ON && Solverptr->itrn_time_now>=Parptr->killsim_time)//-kill命令运行时会开启killsim==on 并且设置killsim_time
    // if not (autosteady state is requested and we're steady)
    && !(Statesptr->steadycheck==ON && Solverptr->t>=Parptr->steadyTotal && !(steadyCount >= 10))//-steady命令运行时会开启steadycheck==on
    )
  {

    // possibly move these variables to global
    iterateq_step();

  }
  //END main ITERATIONS


  return;
}

//----------------------------------------------------------------------------
// SUM Qs INTO A CELL AND UPDATE DEPTH ACCORDINGLY
void UpdateH(States *Statesptr, Pars *Parptr, Solver *Solverptr,BoundCs *BCptr,ChannelSegmentType *ChannelSegments,Arrays *Arrptr)
{
  int i,j,pi;
  double *qxptr0,*qyptr0,*qyptr1,*hptr;
  //qxptr1=qxptr0+1;
  int *mptr;
  double dV,himp,qtmp;
  double dAPorTemp;

    // Insert point sources ((MT) moved before dV check, in case flow out of cell results in negative H reset and the inflow would prevent this)
  // H point sources moved back to after UpdateH (JCN)
  BCptr->Qpoint=0.0;
  for(pi=0;pi<BCptr->numPS;pi++)
  {
    if(BCptr->PS_Ident[pi]==4) // QFIX
    { 
      //H+=QFIX*dx/dA
      Arrptr->H[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz]+=BCptr->PS_Val[pi]*Parptr->dx*Solverptr->Tstep/Parptr->dA;//这一步决定了bci中的流量单位是㎡/s
      //Qpoint+=QFIX*dx
      BCptr->Qpoint+=BCptr->PS_Val[pi]*Parptr->dx;
    }
    if(BCptr->PS_Ident[pi]==5) // QVAR
    { //get QVAR
      //若PS_Ident[i]==5，则PS_Val[pi]记录的不是值，而是相应QVAR序列在BCVarlist中的索引
      qtmp=InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]],Solverptr->t);//重要，InterBC函数根据t来读取QVAR数组中的值
      //H+=QVAR*dx/dA
      Arrptr->H[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz]+=qtmp*Parptr->dx*Solverptr->Tstep/Parptr->dA;
      //Qpoint+=QVAR*dx
      BCptr->Qpoint+=qtmp*Parptr->dx;
    }
  }

  // Calculate dV (+ve => outflow) and update NewH
#pragma omp parallel for private( i,qxptr0,qyptr0,qyptr1,hptr,mptr,dV,dAPorTemp) num_threads(Solverptr->ThreadNum)
  for(j=0;j<Parptr->ysz;j++)
  {
    qxptr0=Arrptr->Qx+j*(Parptr->xsz+1);
    qyptr0=Arrptr->Qy+j*(Parptr->xsz+1);
    qyptr1=Arrptr->Qy+(j+1)*(Parptr->xsz+1);
    hptr=Arrptr->H+j*Parptr->xsz;
    mptr=Arrptr->ChanMask+j*Parptr->xsz;
    for(i=0;i<Parptr->xsz;i++)
    {
      if(*mptr==-1)//判断是否是河流 -1应该代表不是
      {
        if(Statesptr->porosity==ON)
        {
          dV=Solverptr->Tstep*(*qxptr0-*(qxptr0+1)+*qyptr0-*qyptr1);
          dAPorTemp=PorArea(i,j,Parptr,Arrptr);
          if(dAPorTemp==0.0) (*hptr)+=0.0;
          else (*hptr)+=dV/dAPorTemp;

          if(*hptr<0) *hptr=0;
        }
        else
        {
          dV=Solverptr->Tstep*(*qxptr0-*(qxptr0+1)+*qyptr0-*qyptr1);//dV=(Qwest+Qnorth-Qeast-Qsouth)*Tstep
          (*hptr)+=dV/Parptr->dA;
          if(*hptr<0) *hptr=0;
        }
      }
      qxptr0++;
      qyptr0++;
      qyptr1++;
      hptr++;
      mptr++;
    }
  }

  // Point source HVAR and HFIX
  for(pi=0;pi<BCptr->numPS;pi++)
  {
    if(BCptr->PS_Ident[pi]==2) // HFIX
    {
      himp=BCptr->PS_Val[pi]-Arrptr->DEM[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz];
      if(himp<0.0) himp=0.0;
      BCptr->Qpoint+=(himp-Arrptr->H[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz])*Parptr->dA/Solverptr->Tstep;
      Arrptr->H[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz]=himp;
    }
    if(BCptr->PS_Ident[pi]==3) // HVAR
    {
      himp=InterpBC(BCptr->BCVarlist[(int)BCptr->PS_Val[pi]],Solverptr->t)-Arrptr->DEM[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz];
      if(himp<0.0) himp=0.0;
      BCptr->Qpoint+=(himp-Arrptr->H[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz])*Parptr->dA/Solverptr->Tstep;
      Arrptr->H[BCptr->xpi[pi]+BCptr->ypi[pi]*Parptr->xsz]=himp;
    }
  }

  return;
}
