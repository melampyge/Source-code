/*		model.h			*/
//////////////////////////////////////////
// model specific function definitions	//
//////////////////////////////////////////

#ifndef __MODEL_H_
#define __MODEL_H_

#include <cstring>
#include <cmath>

#include "parameter.h"
#include "basics.h"
#include "bc_func.h"
#include "debug.h"
#include "localstress.h"


#ifdef __DEBUG
#include <iomanip>
#include <iostream>
#endif



void calcforces(long, long, long, int);
void calcintracellforces(long, int);

void calcfd(long, long, int);
void calcfd(long, long, double, double, double, double, double, int);
void calcintracellfd(long, int);
void calcintracellfd(long, double, double, double, double, double, int);
void calcrandf(long, long, long);




void calcbackgroundforces(long);

void calcbackgroundfd(long);



////////////////////////////////////
// inline function definitions	  //
////////////////////////////////////


////////////////////////////////////
// double dist{x,y,z}(...)
// ---------------
// calculate {x,y,z} distance between i and j
// ===
// Parameters:
//   unsigned long i: index  of cell particle 1
//   unsigned long j: index  of cell particle 2
// ===
// Returns:
//   double: distance in {x,y,z} direction between particles
inline double distx(long i, long j)
{
  //////////////////////////////////
  // calculate distance between i and j
#ifdef PBC
  return (_rx[j] - _rx[i] - LX*_nearbyint((_rx[j] - _rx[i])/LX));
#endif
#if defined (BBC) || defined (RBC)
  return (_rx[j] - _rx[i]);
#endif
#ifdef HBC
  return (_rx[j] - _rx[i] - LX*_nearbyint((_rx[j] - _rx[i])/LX));
#endif
}

inline double disty(long i, long j)
{
  //////////////////////////////////
  // calculate distance between i and j
#ifdef PBC
  return (_ry[j] - _ry[i] - LY*_nearbyint((_ry[j] - _ry[i])/LY));
#endif
#if defined (BBC) || defined (RBC)
  return (_ry[j] - _ry[i]);
#endif
#ifdef HBC
  return (_ry[j] - _ry[i] - LY*_nearbyint((_ry[j] - _ry[i])/LY));
#endif
}

inline double distz(long i, long j)
{
  //////////////////////////////////
  // calculate distance between i and j
#ifdef PBC
  return (_rz[j] - _rz[i] - LZ*_nearbyint((_rz[j] - _rz[i])/LZ));
#endif
#if defined (BBC) || defined (RBC)
  return (_rz[j] - _rz[i]);
#endif
#ifdef HBC
  return (_rz[j] - _rz[i]);
#endif
}

inline double zpos(long i)
{
  return _rz[i];
}

////////////////////////////////////
// double distsq(...)
// ---------------
// calculate distance squared between i and j
// ===
// Parameters:
//   unsigned long i: index of cell particle 1
//   unsigned long j: index  of cell particle 2
// ===
// Returns:
//   double: distance between particles squared
inline double distsq(long i, long j)
{
  return distx(i,j)*distx(i,j)+disty(i,j)*disty(i,j)+distz(i,j)*distz(i,j);
}



////////////////////////////////////
// double fold{x,y,z}(...)
// ---------------
// fold {x,y,z} coordinates
// ===
// Parameters:
//   unsigned long i: index  of cell particle 1
//   unsigned long j: index  of cell particle 2
// ===
// Returns:
//   double: distance in {x,y,z} direction between particles
inline double foldx(double x)
{
#ifdef PBC
  return (x - LX*_floor(x/LX));
#endif
#if defined (BBC) || defined (RBC)
  return (x);
#endif
#ifdef HBC
  return (x - LX*_floor(x/LX));
#endif
}

inline double foldy(double y)
{
#ifdef PBC
  return (y - LY*_floor(y/LY));
#endif
#if defined (BBC) || defined (RBC)
  return (y);
#endif
#ifdef HBC
  return (y - LY*_floor(y/LY));
#endif
}

inline double foldz(double z)
{
#ifdef PBC
  return (z - LZ*_floor(z/LZ));
#endif
#if defined (BBC) || defined (RBC)
  return (z);
#endif
#ifdef HBC
  return (z);
#endif
}


#endif
