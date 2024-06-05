/* Copyright(c) 2017, The Regents of the University of California, Davis.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met :
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and / or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <algorithm>
#include <assert.h>
#include <cfloat>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>


#ifndef HAVE_VERT
#define HAVE_VERT
#define MAX_CONNECT 20 // the max valence for a triangle
struct vert{
    double x[3];//x,y,z
    int connect[MAX_CONNECT];
    double dih_ang;
};
#endif

struct tri{
    int id[3];
    int neighbour[3];
};

#include "KdTree.h"  // TODO - try nanoflann.hpp here
#include "ThreeDRotation.h"

struct BoundBox{
    double xmin, xmax, ymin, ymax, zmin, zmax, lx, ly, lz;
};

class CMeshImp
{
    public:
	CMeshImp(int numVert, double **Verts, int numTri, int**Tris);
	~CMeshImp();


	//** Non-obtuse Remeshing **//
	void NonobtuseRemeshing(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose);
	void NonobtuseRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose);
	void AcuteRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose);
	void SmallAngleElimination_InterleaveOpt(double targetMinAngle, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double maxAngleAllow, bool verbose);


	//** Sifting/Simplification **//
	void Simp(int targetNumSamples, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose);


	//** Postprocessing **//
	void WriteStatToFile(std::string, int, vert*, bool obt);
	void DisplayStat(int, vert*, bool obt);
	void GetMeshOBJ(std::string filename, bool obtuse, std::string filename_obt);


	void PerturbVertices();

    private:

	bool isObtuse;
	bool isAcute;

	//Classes
	KdTree surfaceKdTree; // the input surface kd tree
	BoundBox myBoundingBox;

	//Functions
	void ErrWarnMessage(size_t, std::string, size_t);
	void FindNeighbourTriangles(int**&);
	bool ObtuseHead(int, int&, int&);
	bool AcuteHead(int ip, int&ip1, int&ip2, double measureAngle);
	bool TooSmallHeadAngle(int ip, int&ip1, int&ip2, double measureAngle);
	void InitialSortAroundVertex(int**);
	void ScaleInputMesh();
	void GetFanTriangle();

	//Variables
	size_t MaxNumVert;
	int numVert_org, numTri_org, numVert_imp;
	vert *Vert_org;//the original surface
	tri *Tri_org;//the original surface triangles
	vert *Vert_imp;//the improved/modified surface
	double scale_factor;
	std::vector<int>indices;

};

#ifdef CMeshImpl

#define MAXPOOL 2000


#define PI 3.14159265358979323846264338327950288
#define TwoPI 6.2831853071795862

#define _tol 10E-8 //this tolerance is based on the unit box. may change for different configuration
#define _tol_sq _tol*_tol
#define _tol_sq_circ (_tol+1.0)*(_tol+1.0) //use this when compare distance^2 with radius^2
					   //basically to check if the a point is inside a circle

#define RadToDeg 57.295779513078550 //180.0/PI (muliply this by the radian angle to convert to degree)

#define _ang_tol_2 5.0 //angle tolerance
#define _ang_tol 3.0 //angle tolerance
#define SIZE_T_MAX ((size_t) -1)

//Structs

struct sphere
{
    double x[4];//x,y,z and r
};
struct plane
{
    double x[4]; //a, b,c and d in the scalar equation of a plane

};

struct stats
{
    double min_ang, max_ang, min_dih_ang, max_dih_ang, min_tri_quality, max_tri_quality;
};

template<typename T>
int GetIndex(T entry, T*list)
{
    //get the index of entry in array (list) of type T
    //if entry is not there, return -1
    //list store the total number of its entries in first element [0]
    for (int i = 1; i <= list[0]; i++){
	if (list[i] == entry){
	    return i;
	}
    }
    return -1;
}

template <typename T>
bool FindCommonElements(T*arr1, T*arr2, T commonEle[])
{
    //store the common elements in arr1 and arr2 in commonEle
    //return true if it is more than one elements

    commonEle[0] = 0;
    int num_arr1(arr1[0]), num_arr2(arr2[0]);
    for (int i = 1; i <= num_arr1; i++){
	for (int j = 1; j <= num_arr2; j++){
	    if (arr1[i] == arr2[j]){
		commonEle[++commonEle[0]] = arr1[i];
	    }
	}
    }

    if (commonEle[0] > 0){ return true; }
    else{ return false; }
}

template <typename T>
bool FindCommonElements_SkipList(T*arr1, T*arr2, T*skip_list, T commonEle[])
{
    //store the common elements in arr1 and arr2 in commonEle
    //return true if it is more than one elements

    commonEle[0] = 0;
    int num_arr1(arr1[0]), num_arr2(arr2[0]), num_skip(skip_list[0]);

    for (int i = 1; i <= num_arr1; i++){
	if (GetIndex(arr1[i], skip_list) >= 0){ continue; }
	for (int j = 1; j <= num_arr2; j++){
	    if (arr1[i] == arr2[j]){
		commonEle[++commonEle[0]] = arr1[i];
	    }
	}
    }

    if (commonEle[0] > 0){ return true; }
    else{ return false; }
}

template <typename T>
T FindCommonElement_SkipList(T*arr1, T*arr2, T*skip_list)
{

    //return the single common elements bteween arr1 and arr2
    //make no check if there is more than one
    int num_arr1(arr1[0]), num_arr2(arr2[0]), num_skip(skip_list[0]);

    for (int i = 1; i <= num_arr1; i++){
	if (GetIndex(arr1[i], skip_list) >= 0){ continue; }
	for (int j = 1; j <= num_arr2; j++){
	    if (arr1[i] == arr2[j]){
		return arr1[i];
	    }
	}
    }

    return -1;
}

template <typename T>
T FindCommonElement_SkipList(T*arr1, T*arr2, T skip_element)
{

    //return the single common elements bteween arr1 and arr2
    //make no check if there is more than one
    int num_arr1(arr1[0]), num_arr2(arr2[0]);

    for (int i = 1; i <= num_arr1; i++){
	if (arr1[i] == skip_element){ continue; }
	for (int j = 1; j <= num_arr2; j++){
	    if (arr1[i] == arr2[j]){
		return arr1[i];
	    }
	}
    }

    return -1;
}

template <typename T>
inline bool FindDuplication(T*arr1)
{
    //find if there is any duplication in arr1
    //arr1 lenght is stored in first element of it
    if (arr1[0] < 2){ return false; }
    for (int i = 1; i <= arr1[0] - 1; i++){
	for (int j = i + 1; j <= arr1[0]; j++){
	    if (arr1[i] == arr1[j]){ return true; }
	}
    }
    return false;
}

template <typename T>
inline void RemoveNodeFromList(T*list, T entry)
{
    //remove entry from list while preseving its sorted format
    int tmp, tmp1, d;
    bool find(false);
    for (d = 1; d <= list[0]; d++){
	if (list[d] == entry){ find = true; break; }
    }
    //entry is not in the list
    if (!find){ return; }
    d = list[0];
    tmp = list[d];
    while (tmp != entry){
	tmp1 = list[d - 1];
	list[d - 1] = tmp;
	tmp = tmp1;
	d--;
    }
    list[0]--;
}

template <typename T>
void AddEntrySortedList(T entry, T* list, T prv, T nxt)
{
    //add entry to the sorted list
    //such that the new entry fits between prv_entry and nxt_entry
    //length of list is stored in list[0]

    if ((prv == list[list[0]] && nxt == list[1]) || (nxt == list[list[0]] && prv == list[1])){
	//it fits to be in the last or the start of list
	list[++list[0]] = entry;
	return;
    }

    if (GetIndex(prv, list) > GetIndex(nxt, list)){
	std::swap(prv, nxt);
    }

    int i(++list[0]);
    while (list[i - 1] != prv){
	list[i] = list[i - 1];
	--i;
    }
    list[i] = entry;
}

inline void NormalizeVector(double&vector_x, double&vector_y, double&vector_z)
{
    double nn = sqrt(vector_x*vector_x + vector_y*vector_y + vector_z*vector_z);
    vector_x /= nn; vector_y /= nn; vector_z /= nn;
}

inline void Cross(double xv1, double yv1, double zv1, double xv2, double yv2, double zv2, double&xx, double&yy, double&zz)
{
    xx = yv1*zv2 - zv1*yv2;
    yy = zv1*xv2 - xv1*zv2;
    zz = xv1*yv2 - yv1*xv2;
}
inline double Dot(double xv1, double yv1, double zv1, double xv2, double yv2, double zv2)
{
    return xv1*xv2 + yv1*yv2 + zv1*zv2;
}
inline void PlaneNorm(double ux, double uy, double uz, double vx, double vy, double vz, double&x, double&y, double&z)
{
    //normal vector to plane defined by two vectors

    x = uy*vz - uz*vy;
    y = uz*vx - ux*vz;
    z = ux*vy - uy*vx;

    double len = sqrt(x*x + y*y + z*z);
    x /= len;
    y /= len;
    z /= len;
}

inline double Angle360Vectors(double dxba, double dyba, double dzba, double dxca, double dyca, double dzca, double norm_x, double norm_y, double norm_z)
{
    // return angle in degrees between b.a.c from (0,2PI)
    //(norm_x,norm_y,norm_z) are the perpendicular normal vector to the triangle (bac)

    double dot, pcross, theta;

    double i, j, k; // cross(ab,ac)

    dot = dxba*dxca + dyba*dyca + dzba*dzca; // dot(ab,ac)

    Cross(dxba, dyba, dzba, dxca, dyca, dzca, i, j, k);

    pcross = i*norm_x + j*norm_y + k*norm_z;

    theta = atan2(pcross, dot);

    if (theta < 0){ theta += TwoPI; }

    return theta *RadToDeg;//convert to deg


}

inline double AngleVectVect(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double dot = x1*x2 + y1*y2 + z1*z2;

    if (dot == 0.0) { return 1.5707963268; }//90 deg

    double angle = dot / sqrt((x1*x1 + y1*y1 + z1*z1) * (x2*x2 + y2*y2 + z2*z2));

    if (angle>1.0){ return 0.0; }
    if (angle<-1.0){ return PI; }

    return acos(angle);
}

inline double Dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double dx, dy, dz;
    dx = x1 - x2;
    dy = y1 - y2;
    dz = z1 - z2;
    dx *= dx;
    dy *= dy;
    dz *= dz;

    return dx + dy + dz;

}

inline double TriCircumcenter3d(double xa, double ya, double za,
	double xb, double yb, double zb,
	double xc, double yc, double zc,
	double&x_cir, double&y_cir, double&z_cir)
    /*(double*a,double*b,double*c,double*circumcenter,double*xi,double*eta)*/
{
    /*****************************************************************************/
    /*
    /*  tricircumcenter3d()   Find the circumcenter of a triangle in 3D.         */
    /*                                                                           */
    /*  The result is returned both in terms of xyz coordinates and xi-eta       */
    /*  coordinates, relative to the triangle's point `a' (that is, `a' is       */
    /*  the origin of both coordinate systems).  Hence, the xyz coordinates      */
    /*  returned are NOT absolute; one must add the coordinates of `a' to        */
    /*  find the absolute coordinates of the circumcircle.  However, this means                                              */
    /*  that the result is frequently more accurate than would be possible if    */
    /*  absolute coordinates were returned, due to limited floating-point        */
    /*  precision.  In general, the circumradius can be computed much more       */
    /*  accurately.                                                              */
    /*                                                                           */
    /*  The xi-eta coordinate system is defined in terms of the triangle.        */
    /*  Point `a' is the origin of the coordinate system.  The edge `ab' extends */
    /*  one unit along the xi axis.  The edge `ac' extends one unit along the    */
    /*  eta axis.  These coordinate values are useful for linear interpolation.  */
    /*                                                                           */
    /*  If `xi' is NULL on input, the xi-eta coordinates will not be computed.   */
    /*                                                                           */
    /*****************************************************************************/
    //http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    //http://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
    double xba, yba, zba, xca, yca, zca;
    double balength, calength;
    double xcrossbc, ycrossbc, zcrossbc;
    double denominator;
    double xcirca, ycirca, zcirca;

    /* Use coordinates relative to point `a' of the triangle. */
    xba = xb - xa;
    yba = yb - ya;
    zba = zb - za;
    xca = xc - xa;
    yca = yc - ya;
    zca = zc - za;
    /* Squares of lengths of the edges incident to `a'. */
    balength = xba * xba + yba * yba + zba * zba;
    calength = xca * xca + yca * yca + zca * zca;

    //	/* Cross product of these edges. */
    //#ifdef EXACT
    //	/* Use orient2d() from http://www.cs.cmu.edu/~quake/robust.html     */
    //	/*   to ensure a correctly signed (and reasonably accurate) result, */
    //	/*   avoiding any possibility of division by zero.                  */
    //	xcrossbc = orient2d(b[1], b[2], c[1], c[2], a[1], a[2]);
    //	ycrossbc = orient2d(b[2], b[0], c[2], c[0], a[2], a[0]);
    //	zcrossbc = orient2d(b[0], b[1], c[0], c[1], a[0], a[1]);
    //#else
    //	/* Take your chances with floating-point roundoff. */
    //	xcrossbc = yba * zca - yca * zba;
    //	ycrossbc = zba * xca - zca * xba;
    //	zcrossbc = xba * yca - xca * yba;
    //#endif

    /* Take your chances with floating-point roundoff. */
    xcrossbc = yba * zca - yca * zba;
    ycrossbc = zba * xca - zca * xba;
    zcrossbc = xba * yca - xca * yba;

    /* Calculate the denominator of the formulae. */
    denominator = 0.5 / (xcrossbc * xcrossbc + ycrossbc * ycrossbc + zcrossbc * zcrossbc);

    /* Calculate offset (from `a') of circumcenter. */
    xcirca = ((balength * yca - calength * yba) * zcrossbc -
	    (balength * zca - calength * zba) * ycrossbc) * denominator;
    ycirca = ((balength * zca - calength * zba) * xcrossbc -
	    (balength * xca - calength * xba) * zcrossbc) * denominator;
    zcirca = ((balength * xca - calength * xba) * ycrossbc -
	    (balength * yca - calength * yba) * xcrossbc) * denominator;
    x_cir = xcirca + xa;
    y_cir = ycirca + ya;
    z_cir = zcirca + za;

    //if (xi != (double *) NULL) {
    /* To interpolate a linear function at the circumcenter, define a     */
    /*   coordinate system with a xi-axis directed from `a' to `b' and    */
    /*   an eta-axis directed from `a' to `c'.  The values for xi and eta */
    /*   are computed by Cramer's Rule for solving systems of linear      */
    /*   equations.                                                       */

    /* There are three ways to do this calculation - using xcrossbc, */
    /*   ycrossbc, or zcrossbc.  Choose whichever has the largest    */
    /*   magnitude, to improve stability and avoid division by zero. */
    /*if (((xcrossbc >= ycrossbc) ^ (-xcrossbc > ycrossbc)) &&
      ((xcrossbc >= zcrossbc) ^ (-xcrossbc > zcrossbc))) {
     *xi = (ycirca * zca - zcirca * yca) / xcrossbc;
     *eta = (zcirca * yba - ycirca * zba) / xcrossbc;
     } else if ((ycrossbc >= zcrossbc) ^ (-ycrossbc > zcrossbc)) {
     *xi = (zcirca * xca - xcirca * zca) / ycrossbc;
     *eta = (xcirca * zba - zcirca * xba) / ycrossbc;
     } else {
     *xi = (xcirca * yca - ycirca * xca) / zcrossbc;
     *eta = (ycirca * xba - xcirca * yba) / zcrossbc;
     }
     }*/

    double len1, dx, dy, dz;
    dx = xa - x_cir;
    dy = ya - y_cir;
    dz = za - z_cir;
    len1 = dx*dx + dy*dy + dz*dz;
    if (len1>2.0){ return len1; }

    /*#ifdef  debug

      dx = xb - x_cir;
      dy = yb - y_cir;
      dz = zb - z_cir;
      len2 = dx*dx + dy*dy + dz*dz;

      dx = xc - x_cir;
      dy = yc - y_cir;
      dz = zc - z_cir;
      len3 = dx*dx + dy*dy + dz*dz;

      if (abs(len1 - len2)>_tol || abs(len3 - len2)>_tol || abs(len1 - len3)>_tol){
    //return 1000;
    cout << "Error at TriCircumcenter3d()..!! " << endl;
    //system("pause");
    }
#endif*/

    return len1;


}

inline double TriTriNormalAngle3(double xip, double yip, double zip, double*ip1, double*ip2, double*ip3)
{
    //return the angles between the normal of the planes containing the two triangles sharing an edge ip-ip1
    //watch out for how you call it. the order of the points is curcial for return the correct angle

    double x_n1, y_n1, z_n1, x_n2, y_n2, z_n2, x_np, y_np, z_np, angle360;

    Cross(ip1[0] - xip, ip1[1] - yip, ip1[2] - zip,
	    ip2[0] - xip, ip2[1] - yip, ip2[2] - zip,
	    x_n1, y_n1, z_n1);
    Cross(ip3[0] - xip, ip3[1] - yip, ip3[2] - zip,
	    ip1[0] - xip, ip1[1] - yip, ip1[2] - zip,
	    x_n2, y_n2, z_n2);

    PlaneNorm(x_n1, y_n1, z_n1, x_n2, y_n2, z_n2, x_np, y_np, z_np);

    angle360 = Angle360Vectors(x_n1, y_n1, z_n1, x_n2, y_n2, z_n2, x_np, y_np, z_np);

    return angle360;
}

inline double TriTriNormalAngle(double*ip, double*ip1, double*ip2, double*ip3)
{
    return TriTriNormalAngle3(ip[0], ip[1], ip[2], ip1, ip2, ip3);
}

inline double PointTriangleDistance(double xp1, double yp1, double zp1,// triangle head1
	double xp2, double yp2, double zp2,// triangle head2
	double xp3, double yp3, double zp3,// triangle head3
	double xp, double yp, double zp, // my point
	double&x_new, double&y_new, double&z_new)// my projection
{
    // http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
    // http://www.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d/content/pointTriangleDistance.m
    //        ^t
    //  \     |
    //   \reg2|
    //    \   |
    //     \  |
    //      \ |
    //       \|
    //        *P2
    //        |\
    //        | \
    //  reg3  |  \ reg1
    //        |   \
    //        |reg0\
    //        |     \
    //        |      \ P1
    // -------*-------*------->s
    //        |P0      \
    //  reg4  | reg5    \ reg6

    /*size_t p1(_surface_tri2[tri][0]),p2(_surface_tri2[tri][1]),p3(_surface_tri2[tri][2]);
      double xp1(_surface[p1][0]),yp1(_surface[p1][1]),zp1(_surface[p1][2]),
      xp2(_surface[p2][0]),yp2(_surface[p2][1]),zp2(_surface[p2][2]),
      xp3(_surface[p3][0]),yp3(_surface[p3][1]),zp3(_surface[p3][2]);*/

    double xE0(xp2 - xp1), yE0(yp2 - yp1), zE0(zp2 - zp1),
	   xE1(xp3 - xp1), yE1(yp3 - yp1), zE1(zp3 - zp1),
	   xD(xp1 - xp), yD(yp1 - yp), zD(zp1 - zp);

    double a = Dot(xE0, yE0, zE0, xE0, yE0, zE0);
    double b = Dot(xE0, yE0, zE0, xE1, yE1, zE1);
    double c = Dot(xE1, yE1, zE1, xE1, yE1, zE1);
    double d = Dot(xE0, yE0, zE0, xD, yD, zD);
    double e = Dot(xE1, yE1, zE1, xD, yD, zD);
    double f = Dot(xD, yD, zD, xD, yD, zD);

    double det = a*c - b*b;
    double s = b*e - c*d;
    double t = b*d - a*e;
    double invDet, dist, numer, denom, tmp0, tmp1;
    //Terible tree of conditionals to determine in which region of the diagram
    // shown above the projection of the point into the triangle-plane lies.

    if (s + t <= det){
	if (s<0){
	    if (t<0){
		//region 4
		if (d < 0){
		    t = 0;
		    if (-1.0*d >= a){
			s = 1;
			dist = a + 2.0*d + f;
		    }
		    else{
			s = -1.0*d / a;
			dist = d*s + f;
		    }
		}
		else{
		    s = 0;
		    if (e >= 0){
			t = 0;
			dist = f;
		    }
		    else{
			if (-e >= c){
			    t = 1;
			    dist = c + 2.0*e + f;
			}
			else{
			    t = -1.0*e / c;
			    dist = e*t + f;
			}

		    }

		}

	    }
	    else
	    {
		//region 3
		s = 0;
		if (e >= 0){
		    t = 0;
		    dist = f;
		}
		else{
		    if (-1.0*e >= c){
			t = 1;
			dist = c + 2.0*e + f;
		    }
		    else{
			t = -1.0*e / c;
			dist = e*t + f;
		    }
		}

	    }
	}
	else if (t<0){
	    //region 5
	    t = 0;
	    if (d >= 0){
		s = 0;
		dist = f;
	    }
	    else{
		if (-1.0*d >= a){
		    s = 1;
		    dist = a + 2.0*d + f;
		}
		else{
		    s = -1.0*d / a;
		    dist = d*s + f;
		}
	    }

	}
	else{
	    //region 0
	    invDet = 1.0 / det;
	    s = s*invDet;
	    t = t*invDet;
	    dist = s*(a*s + b*t + 2.0*d) +
		t*(b*s + c*t + 2.0*e) + f;
	}
    }
    else{
	if (s<0){
	    //region 2
	    tmp0 = b + d;
	    tmp1 = c + e;
	    if (tmp1 > tmp0){ // minimum on edge s+t=1
		numer = tmp1 - tmp0;
		denom = a - 2.0*b + c;
		if (numer >= denom){
		    s = 1;
		    t = 0;
		    dist = a + 2.0*d + f;
		}
		else{
		    s = numer / denom;
		    t = 1.0 - s;
		    dist = s*(a*s + b*t + 2.0*d)
			+ t*(b*s + c*t + 2 * e) + f;
		}
	    }
	    else{
		s = 0;
		if (tmp1 <= 0){
		    t = 1;
		    dist = c + 2.0*e + f;
		}
		else{
		    if (e >= 0){
			t = 0;
			dist = f;
		    }
		    else{
			t = -1.0*e / c;
			dist = e*t + f;
		    }
		}
	    }
	}
	else if (t<0){
	    //region 6

	    tmp0 = b + e;
	    tmp1 = a + d;
	    if (tmp1 > tmp0){ // minimum on edge s+t=1
		numer = tmp1 - tmp0;
		denom = a - 2.0*b + c;
		if (numer >= denom){
		    t = 1;
		    s = 0;
		    dist = c + 2.0*e + f;
		}
		else{
		    t = numer / denom;
		    s = 1.0 - t;
		    dist = s*(a*s + b*t + 2.0*d)
			+ t*(b*s + c*t + 2.0*e) + f;
		}
	    }
	    else{
		t = 0;
		if (tmp1 <= 0){
		    s = 1;
		    dist = a + 2.0*d + f;
		}
		else{
		    if (d >= 0){
			s = 0;
			dist = f;
		    }
		    else{
			s = -d / a;
			dist = d*s + f;
		    }
		}
	    }
	}
	else{
	    //region 1
	    numer = c + e - b - d;
	    if (numer <= 0){
		s = 0;
		t = 1;
		dist = c + 2.0*e + f;
	    }
	    else{
		denom = a - 2.0*b + c;
		if (numer >= denom){
		    s = 1;
		    t = 0;
		    dist = a + 2.0*d + f;
		}
		else{
		    s = numer / denom;
		    t = 1 - s;
		    dist = s*(a*s + b*t + 2.0*d) +
			t*(b*s + c*t + 2.0*e) + f;
		}
	    }
	}
    }

    x_new = xp1 + s*xE0 + t*xE1;
    y_new = yp1 + s*yE0 + t*yE1;
    z_new = zp1 + s*zE0 + t*zE1;
    if (dist<0){ dist = 0; }
    return dist;

}

class RndNum
{
    public:
	RndNum();
	~RndNum();

	void InitiateRndNum(unsigned long x);
	double RandNumGenerator();


    private:
	int indx;
	double Q[1220];
};

//TODO get different stats at single vertex
//TODO report valence
class Statistics
{
    public:
	Statistics();
	~Statistics();

	//First call this to calc all statistics
	void GetAllStats(int numVert, vert *Vert);

	//Quary functions for the calculated statistics
	double GetMinAngle(){ return min_angle; }
	double GetMaxAngle(){ return max_angle; }
	void GetMinMaxEdgeLen(double&myMinEdgeLen, double&myMaxEdgeLen){ myMinEdgeLen = min_edge_len; myMaxEdgeLen = max_edge_len; }
	void GetMinMaxTriQuality(double&myMinQuality, double&myMaxQuality){ myMinQuality = min_quality; myMaxQuality = max_quality; }
	int GetNumTri(){ return num_tri; }
	int GetNumNonObtuse(){ return num_tri_obtuse; }
	int CalcNumAcute(int numVert, vert *Vert, double measureAngle);
	void DisplayStatistics(FILE *fp, bool obt);

    private:
	double min_dih_ang, max_dih_ang,
	       min_quality, max_quality,
	       min_angle, max_angle,
	       min_edge_len, max_edge_len,
	       average_div_from_60; //average deviation of interior angles from 60 degree

	int num_tri, num_tri_obtuse, num_vertices;


	double getCurv(int, int, vert *);
	void CalcDihAngle(int, vert *);
	void CalcAnglesEdgeLen(int, vert *);
	void CalcTriangelQuality(int, vert *);
	void CalcNumTriangles(int, vert *);
	void CalcNumNonObtuse(int, vert *);
	double SingleTriangleQuality(int, int, int, vert *);
	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
};

class Constraints
{

    public:
	Constraints();
	~Constraints();
	void Reset(vert*Verts, int*nList);
	bool Delaunay(int*nList, int*skipList, vert*Verts, bool loadSkipList, bool isConnected);
	bool DelaunayForcedNotConnected(int*nList, int removed, int*skipList, vert*Verts);
	void MinMaxAngle(int*nList, vert*Verts, double min_ang, double max_ang, double*void_vertex, bool nonObtuse, bool isConnected);
	void MinEdgeLength(int*nList, vert*Verts, double (*sizingfunc)(double xx, double yy, double zz));
	void MinEdgeLength(int*nList, vert*Verts, double r_min_2);
	void Smoothness(int*nList, int*skipList, vert*Verts, double*void_vertex, double dev);
	void Direction(int*ip, int*nList, vert*Verts);
	void Maximality();
	void AttractorInSphere(vert*Verts, int ip, double*void_vertex);
	void RepellerExSphere(vert*Verts, int ip, double*void_vertex);

	bool InsideFeasibleRegion_Vertex(double xx, double yy, double zz);
	bool InsideFeasibleRegion_Triangle(double*tri);
	bool OverlappingInSpheres();
	void SinglePlane(double xn, double yn, double zn, double px, double py, double pz);

	bool *isObtuse = NULL;
	bool *isAcute = NULL;

	//debug
	void DrawExSpheres();
	void DrawInSpheres();


    private:

	int numExSphere, numInSphere, numPlane;
	int numExSphere_size, numInSphere_size, numPlane_size;
	int aux_list[20];

	sphere *ExSphere, *InSphere;
	plane *Plane;
	BoundBox myBox; //specify the bounding box using the nList
			//new vertex should be inside the bounding box
	void SetBoundingBox(vert*Verts, int*nList);
	bool DelaunayNotConnected(int*nList, int*skipList, vert*Verts, bool loadSkipList);
	void ExpandSpheres(int&currentSize, int currentNumSphere, sphere*&mySpheres);
	void ExpandPlanes(int&currentSize, int currentNumPlane, plane*&myPlanes);
	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
	bool IsEmptySphere(double xc, double yc, double zc, double rc_2,
		int*list,
		int*skip_list,
		vert*Verts);
	bool IsInsideBoundingBox(double xx, double yy, double zz);
	//debug
	void SpheresDrawing(std::string filename, int num, sphere *Sphere);
};

class Operator
{
    public:
	//the constraints are set once by the executer
	struct AppConstraints{
	    bool isDelaunay;
	    bool isMinAngle; double MinAngle;
	    bool isMaxAngle; double MaxAngle;
	    bool isNonobtuse;
	    bool isSmooth; double dev;
	    bool isEdgeLen; double MinEdgeLength_sq;
	    //need to look into this in case of function is passed
	    bool isMaximal;
	} constraints;

	Operator(vert*Vert_org, bool *isO, bool *isA);
	~Operator();

	void TriValentRemoval(int&numVert, vert*Verts);

	bool Relocation(int ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer);

	bool Ejection(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool att);

	bool Injection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj);

	bool AggressiveInjection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj);



    private:
	Constraints myConstraints;
	RndNum myRandNum;

	int skipList[20];//the skip list stores the vertices that should be considered as if they are deleted
	int skipListNN[20];//the skip list stores the vertices that should be considered as if they are deleted
			   //we use this when we are attracting/repelling vertices
	double newVertex[3];//coordinates of the new vertex (updated from Sampler)
	int mynList[100];//sorted neighbour list (connected chain)
	int mynListNN[100];//sorted neighbour list for being attracted vertex (unconnected chain)
	double void_vertex[3];
	int apex[4];//for each two consecutive vertices in *ip, find their (correct) shared vertex
	int EdgeToBreak[10];//max 3 edges (3*2 + 1)
	double original_neighbout[100][3];//keeps a copy of the original void corners to be updated in case of failure in repeller or attractor
	int temp_arr[100];//used for various reasons when a temp array is needed

	void AttractorRepeller(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool isAttractor);
	void SetAttractorRepeller(vert*Verts);
	void RevertAttractorRepeller(vert*Verts);

	void StartActivePool(int closedtSurfaceID, double*myVert);
	void StartActivePool(int id, int numLayers);

	void GetFanTriangles(int);
	bool IsDuplicated(int ip1, int ip2, int ip3, int t2);
	bool Sampler(vert*Verts, int*nList, int budget,
		int i_af, int i_rep, int i_bf,
		double(*Optimizer)(vert*Verts, int*nList, double*void_ver, double xx, double yy, double zz));

	void RetrieveCoordinates(int lf, int* cell, double&x0, double&y0, double&z0, double&x1, double&y1, double&z1, double&x2, double&y2, double&z2);
	bool CheckNewVertex(vert*Verts, int*nList, double x_new, double  y_new, double  z_new, bool att);
	bool CheckRefinement(vert*Verts, int*nList, double*tri);
	void RandomPointInTri(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double&x, double&y, double&z);
	void Mover(int ip, vert*Verts);
	void EdgeCollapse(int ip_low, int ip_high, int&numVert, vert*Verts);
	void FaceCollapse(int ip_low, int ip_mid, int ip_high, int&numVert, vert*Verts);
	void Inserter(int&numVert, vert*Verts);
	void RemoveVertex(int ip, int&numVert, vert*Verts);
	bool InspectFeasibleRegion(int*nList, vert*Verts);

	void ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id);
	void TriValent(int ip, vert*Verts, int*list, int skip1, int skip2);
	void GetListSkipVertex(int ip, vert*Verts, int skip1, int skip2, int*list);
	void GetEdgeSortedNeighbourList(int ip1, int ip2, vert*Verts, int*list);
	void GetFaceSortedNeighbourList(int ip1, int ip2, int ip3, vert*Verts, int*list);
	void CycleList(int ip, vert*Verts, int iq, int start_id, int end_id, int*list, bool includeLastEle);
	void CycleList_SkipList(int ip, vert*Verts, int*SkipList, int start_id, int end_id, int*list);


	void DrawActivePool(int lf);
	void DrawnList(vert*Verts);

	vert*Vert_org;//input surface

	//resampling grid
	int num_active;
	int**active_pool, **tmp_active_pool, **tri_pool;
	double **_tar;
	int next_layer[MAXPOOL];
};

RndNum::RndNum()
{
}
RndNum::~RndNum()
{
}

void RndNum::InitiateRndNum(unsigned long x)
{
    assert(sizeof(double) >= 8);

    size_t i;
    size_t qlen = indx = sizeof Q / sizeof Q[0];
    for (i = 0; i < qlen; i++)
	Q[i] = 0;
    double c = 0.0, zc = 0.0,	/* current CSWB and SWB `borrow` */
	   zx = 5212886298506819.0 / 9007199254740992.0,	/* SWB seed1 */
	   zy = 2020898595989513.0 / 9007199254740992.0;	/* SWB seed2 */
    int j;
    double s, t;	 /* Choose 32 bits for x, 32 for y */
    unsigned long /*x = 123456789,*/ y = 362436069; /* default seeds * /
    /* Next, seed each Q[i], one bit at a time, */

    if (x == 0){ x = 123456789; }

    for (i = 0; i < qlen; i++)
    { /* using 9th bit from Cong+Xorshift */
	s = 0.0;
	t = 1.0;
	for (j = 0; j < 52; j++)
	{
	    t = 0.5 * t; /* make t=.5/2^j */
	    x = 69069 * x + 123;
	    y ^= (y << 13);
	    y ^= (y >> 17);
	    y ^= (y << 5);
	    if (((x + y) >> 23) & 1)
		s = s + t; /* change bit of s, maybe */
	}	 /* end j loop */
	Q[i] = s;
    } /* end i seed loop, Now generate 10^9 RandNumGenerator's: */
}
double RndNum::RandNumGenerator()
{
    double c = 0.0, zc = 0.0,	/* current CSWB and SWB `borrow` */
	   zx = 5212886298506819.0 / 9007199254740992.0,	/* SWB seed1 */
	   zy = 2020898595989513.0 / 9007199254740992.0;	/* SWB seed2 */

    /* Takes 14 nanosecs, Intel Q6600,2.40GHz */
    int i, j;
    double t; /* t: first temp, then next CSWB value */
    /* First get zy as next lag-2 SWB */
    t = zx - zy - zc;
    zx = zy;
    double cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
    if (t < 0)
    {
	zy = t + 1.0;
	zc = cc;
    }
    else
    {
	zy = t;
	zc = 0.0;
    }

    /* Then get t as the next lag-1220 CSWB value */
    if (indx < 1220)
	t = Q[indx++];
    else
    { /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
	for (i = 0; i < 1220; i++)
	{
	    j = (i < 30) ? i + 1190 : i - 30;
	    t = Q[j] - Q[i] + c; /* Get next CSWB element */
	    if (t > 0)
	    {
		t = t - cc;
		c = cc;
	    }
	    else
	    {
		t = t - cc + 1.0;
		c = 0.0;
	    }
	    Q[i] = t;
	}	 /* end i loop */
	indx = 1;
	t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
    } /* end else segment; return t-zy mod 1 */
    return ((t < zy) ? 1.0 + (t - zy) : t - zy);
} /* end RandNumGenerator() */

//TODO report non-delaunay triangle
Statistics::Statistics()
{
    min_dih_ang = DBL_MAX;
    max_dih_ang = DBL_MIN;
    min_quality = DBL_MAX;
    max_quality = DBL_MIN;
    min_angle = DBL_MAX;
    max_angle = DBL_MIN;
    min_edge_len = DBL_MAX;
    max_edge_len = DBL_MIN;
    num_tri = 0;
    num_tri_obtuse = 0;
}

Statistics::~Statistics()
{
}

void Statistics::DisplayStatistics(FILE *fp, bool obt)
{
    fprintf(fp, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(fp, "********************************** Statistics ***********************************");
    fprintf(fp, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    fprintf(fp, "\n #Triangle = %i", num_tri);
    fprintf(fp, "\n #Vertex = %i", num_vertices);
    fprintf(fp, "\n");

    if (obt){
	fprintf(fp, "\n #Obtuse triangle = %i", num_tri_obtuse);
	fprintf(fp, "\n Obtuse triangle = %f %%", 100.0*double(num_tri_obtuse) / double(num_tri));
	fprintf(fp, "\n");
    }

    fprintf(fp, "\n Minimum angle = %f", min_angle);
    fprintf(fp, "\n Maximum angle = %f", max_angle);
    fprintf(fp, "\n Average deviation from 60 deg = %f", average_div_from_60);
    fprintf(fp, "\n");

    fprintf(fp, "\n Minimum edge length = %f", sqrt(min_edge_len));
    fprintf(fp, "\n Maximum edge length = %f", sqrt(max_edge_len));
    fprintf(fp, "\n Minimum to Maximum edge length = %f  %%", 100.0*sqrt(min_edge_len / max_edge_len));
    fprintf(fp, "\n");

    fprintf(fp, "\n Minimum triangle quality = %f", min_quality);
    fprintf(fp, "\n Maximum triangle quality = %f", max_quality);
    fprintf(fp, "\n");

    fprintf(fp, "\n Minimum dihedral angle = %f", min_dih_ang);
    fprintf(fp, "\n Maximum dihedral angle = %f", max_dih_ang);
    fprintf(fp, "\n");
}

void Statistics::GetAllStats(int numVert, vert *Vert)
{

    CalcDihAngle(numVert, Vert);
    CalcAnglesEdgeLen(numVert, Vert);
    CalcTriangelQuality(numVert, Vert);
    CalcNumTriangles(numVert, Vert);
    CalcNumNonObtuse(numVert, Vert);
}

//** Num triangles **//
void Statistics::CalcNumTriangles(int numVert, vert *Vert)
{
    num_tri = 0;
    num_vertices = numVert;
    for (int id0 = 0; id0 < numVert; id0++){
	int id1 = Vert[id0].connect[Vert[id0].connect[0]];
	for (int i = 1; i <= Vert[id0].connect[0]; i++){
	    int id2 = Vert[id0].connect[i];
	    num_tri++;
	}
    }

    if (num_tri % 3 != 0){
	ErrWarnMessage(__LINE__, "Statistics::CalcNumTriangles num_tri%3!=0", 1);
    }

    num_tri /= 3;
}

//** Angle and edge length**//
void Statistics::CalcAnglesEdgeLen(int numVert, vert *Vert)
{
    min_angle = 180.0;
    max_angle = 0.0;
    min_edge_len = DBL_MAX;
    max_edge_len = DBL_MIN;
    num_vertices = numVert;
    average_div_from_60 = 0;
    int num_angles = 0;
    for (int id0 = 0; id0 < numVert; id0++){

	int id1 = Vert[id0].connect[Vert[id0].connect[0]];

	for (int i = 1; i <= Vert[id0].connect[0]; i++){
	    int id2 = Vert[id0].connect[i];

	    double angle = AngleVectVect(Vert[id2].x[0] - Vert[id0].x[0], Vert[id2].x[1] - Vert[id0].x[1], Vert[id2].x[2] - Vert[id0].x[2],
		    Vert[id1].x[0] - Vert[id0].x[0], Vert[id1].x[1] - Vert[id0].x[1], Vert[id1].x[2] - Vert[id0].x[2])*RadToDeg;//180.0 / PI;

	    double len = Dist(Vert[id0].x[0], Vert[id0].x[1], Vert[id0].x[2], Vert[id1].x[0], Vert[id1].x[1], Vert[id1].x[2]);

	    min_angle = std::min(min_angle, angle);
	    max_angle = std::max(max_angle, angle);
	    min_edge_len = std::min(min_edge_len, len);
	    max_edge_len = std::max(max_edge_len, len);

	    average_div_from_60 += abs(angle - 60.0);
	    num_angles++;

	    id1 = id2;
	}
    }

    average_div_from_60 /= double(num_angles);

}

//** Non obtuse triangles**//
void Statistics::CalcNumNonObtuse(int numVert, vert *Vert)
{
    num_tri_obtuse = 0;
    num_vertices = numVert;
    for (int id0 = 0; id0 < numVert; id0++){

	int id1 = Vert[id0].connect[Vert[id0].connect[0]];

	for (int i = 1; i <= Vert[id0].connect[0]; i++){
	    int id2 = Vert[id0].connect[i];

	    double angle = AngleVectVect(Vert[id2].x[0] - Vert[id0].x[0], Vert[id2].x[1] - Vert[id0].x[1], Vert[id2].x[2] - Vert[id0].x[2],
		    Vert[id1].x[0] - Vert[id0].x[0], Vert[id1].x[1] - Vert[id0].x[1], Vert[id1].x[2] - Vert[id0].x[2])*RadToDeg; /* 57.295779513078550 = 180.0 / PI*/
	    if (angle > 90.0 + _tol){ num_tri_obtuse++; }

	    id1 = id2;
	}
    }
}

//** Acute triangles**//
int Statistics::CalcNumAcute(int numVert, vert *Vert, double measureAngle)
{
    int num_acute_angle = 0;
    for (int id0 = 0; id0 < numVert; id0++){

	int id1 = Vert[id0].connect[Vert[id0].connect[0]];

	for (int i = 1; i <= Vert[id0].connect[0]; i++){
	    int id2 = Vert[id0].connect[i];

	    double angle = AngleVectVect(Vert[id2].x[0] - Vert[id0].x[0], Vert[id2].x[1] - Vert[id0].x[1], Vert[id2].x[2] - Vert[id0].x[2],
		    Vert[id1].x[0] - Vert[id0].x[0], Vert[id1].x[1] - Vert[id0].x[1], Vert[id1].x[2] - Vert[id0].x[2])*RadToDeg;
	    if (angle > measureAngle + _tol){ num_acute_angle++; }
	    id1 = id2;
	}
    }
    return num_acute_angle;
}

//** Dih angle**//
double Statistics::getCurv(int ip, int numVert, vert *Vert)
{
    //curv is the supplementary angle of the dihderal angle between two triangles
    //Here for vertex ip, we return the average curv between each two triangle in ip triangle fan

    double curve(0), angle;

    int ip3 = Vert[ip].connect[Vert[ip].connect[0]];

    for (int i = 1; i <= Vert[ip].connect[0]; i++){
	int ip1 = Vert[ip].connect[i];
	int ip2 = (i == Vert[ip].connect[0]) ? Vert[ip].connect[1] : Vert[ip].connect[i + 1];

	angle = TriTriNormalAngle(Vert[ip].x, Vert[ip1].x, Vert[ip2].x, Vert[ip3].x);

	//curve = std::max(curve, angle);
	curve += angle;
	ip3 = ip1;
    }
    return curve / double(Vert[ip].connect[0]);
}
void Statistics::CalcDihAngle(int numVert, vert *Vert)
{
    //Calc and store the dih angle at each vertex
    //the dih angle of a vertex is the average dih angle of the vertex triangle fan

    min_dih_ang = DBL_MAX;
    max_dih_ang = DBL_MIN;
    num_vertices = numVert;

    for (int i = 0; i < numVert; i++){
	Vert[i].dih_ang = 180.0 - getCurv(i, numVert, Vert);
	min_dih_ang = std::min(min_dih_ang, Vert[i].dih_ang);
	max_dih_ang = std::max(max_dih_ang, Vert[i].dih_ang);
    }
}

//** Triangle quality**//
double Statistics::SingleTriangleQuality(int ip, int ip1, int ip2, vert *Vert)
{

    double l1, l2, l3, longest, half_perimeter, area, q;
    l1 = sqrt(Dist(Vert[ip].x[0], Vert[ip].x[1], Vert[ip].x[2], Vert[ip1].x[0], Vert[ip1].x[1], Vert[ip1].x[2]));
    l2 = sqrt(Dist(Vert[ip].x[0], Vert[ip].x[1], Vert[ip].x[2], Vert[ip2].x[0], Vert[ip2].x[1], Vert[ip2].x[2]));
    l3 = sqrt(Dist(Vert[ip2].x[0], Vert[ip2].x[1], Vert[ip2].x[2], Vert[ip1].x[0], Vert[ip1].x[1], Vert[ip1].x[2]));

    longest = std::max(l1, l2);
    longest = std::max(longest, l3);

    half_perimeter = 0.5*(l1 + l2 + l3);

    area = sqrt(half_perimeter*(half_perimeter - l1)*(half_perimeter - l2)*(half_perimeter - l3));
    q = (3.4641016151377548 * area) / (half_perimeter*longest); //3.4641016151377548 = 6.0/sqrt(3.0)
								//q = (6.0*area) / (sqrt(3.0)*half_perimeter*longest);

    return q;


}
void Statistics::CalcTriangelQuality(int numVert, vert *Vert)
{
    min_quality = DBL_MAX;
    max_quality = DBL_MIN;
    num_vertices = numVert;
    for (int ip = 0; ip < numVert; ip++){
	int ip2 = Vert[ip].connect[Vert[ip].connect[0]];
	for (int i = 1; i <= Vert[ip].connect[0]; i++){
	    int ip1 = Vert[ip].connect[i];
	    if (ip < ip1 && ip < ip2){//just do the triangle once
		double qu = SingleTriangleQuality(ip, ip1, ip2, Vert);
		min_quality = std::min(min_quality, qu);
		max_quality = std::max(max_quality, qu);
	    }

	    ip2 = ip1;
	}
    }
}

void Statistics::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{

    if (mess_id == 0){
	std::cerr << "\nError::line(" << lineNum << ")-->>" << message << std::endl;
	system("pause");
    }
    else{
	std::cerr << "\nWarning::line(" << lineNum << ")-->>" << message << std::endl;
	system("pause");
    }
}

Constraints::Constraints()
{
    numPlane = numExSphere = numInSphere = 0;
    numExSphere_size = numInSphere_size = numPlane_size = 999;
    ExSphere = new sphere[1000];
    InSphere = new sphere[1000];
    Plane = new plane[1000];
}
//TODO add maximality constraints as 3d points (on the original input surface) to be covered by the new vertex
//TODO add Hausdorff dist constriants


Constraints::~Constraints(){}

void Constraints::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{

    //mess_id =0 for error (exit)
    //otherwise, it is a warning (pause)

    if (mess_id == 0){
	fprintf(stderr, "\nError::line(%d)-->>%s", lineNum, message.c_str());
	system("pause");
    }
    else{
	fprintf(stderr, "\nWarning::line(%d)-->>%s", lineNum, message.c_str());
	system("pause");
    }
}

void Constraints::Reset(vert*Verts, int*nList)
{
    numPlane = numExSphere = numInSphere = 0;
    SetBoundingBox(Verts, nList);
}

//**** Construction
void Constraints::SetBoundingBox(vert*Verts, int*nList)
{
    //set the boundinng using nList
    //if you don't want the bounding box call this with NULL
    if (nList == NULL){
	myBox.xmin = myBox.ymin = myBox.zmin = DBL_MAX;
	myBox.xmax = myBox.ymax = myBox.zmax = DBL_MIN;
    }
    else{
	myBox.xmin = myBox.ymin = myBox.zmin = DBL_MAX;
	myBox.xmax = myBox.ymax = myBox.zmax = DBL_MIN;
	for (int i = 1; i <= nList[0]; i++){
	    int ip = nList[i];
	    myBox.xmin = std::min(myBox.xmin, Verts[ip].x[0]);
	    myBox.ymin = std::min(myBox.ymin, Verts[ip].x[1]);
	    myBox.zmin = std::min(myBox.zmin, Verts[ip].x[2]);

	    myBox.xmax = std::max(myBox.xmax, Verts[ip].x[0]);
	    myBox.ymax = std::max(myBox.ymax, Verts[ip].x[1]);
	    myBox.zmax = std::max(myBox.zmax, Verts[ip].x[2]);
	}

	myBox.lx = myBox.xmax - myBox.xmin;
	myBox.ly = myBox.ymax - myBox.ymin;
	myBox.lz = myBox.zmax - myBox.zmin;
    }
}
bool Constraints::Delaunay(int*nList, int*skipList, vert*Verts, bool loadSkipList, bool isConnected)
{
    if (!isConnected){
	return DelaunayNotConnected(nList, skipList, Verts, loadSkipList);
    }
    //TODO Delaunay_Extend where the list from which we calc ex sphere and in shperes are decoupled

    //nList is the sort list of neighbours around the a void
    //1) find the exclusion spheres such that a new vertex should lie outside to maintain delaunayness of
    //untouched/surounding triangles
    //2) find the inclusion spheres such that a new vertex should lie within to create delaunay-accurate set
    //of triangles
    //skipList contians current vertices that should be negelected

    //start by loading skipList to aux_list
    //for insert, we always wanna consider the vertices in skipList
    if (loadSkipList){
	for (int i = 1; i <= skipList[0]; i++){
	    aux_list[i + 3] = skipList[i];
	}
	aux_list[0] = skipList[0] + 3;
    }
    else{
	aux_list[0] = 3;
    }


    if (nList[0] == 3){
	//this is easy because this result into one InSphere (that connect these three vertices)
	//and three ExSpheres (no looping needed)

	int ap(nList[1]), n1(nList[2]), n2(nList[3]), ap2;
	InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
	numInSphere++;

	//have to expand skipList
	skipList[++skipList[0]] = ap;
	skipList[++skipList[0]] = n1;
	skipList[++skipList[0]] = n2;

	//first exsphere
	ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n2].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;

	//second exsphere
	ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;

	//third exsphere
	ap2 = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;
	skipList[0] -= 3;
	return  true;
    }




    int n2 = nList[nList[0]];
    for (int i = 1; i <= nList[0]; i++){
	int ap = nList[i];
	int n1 = (i == nList[0]) ? nList[1] : nList[i + 1];

	int ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n2].connect, skipList);
	//find the exclusion sphere (circumsphere of triangle n2-ap-ap2)
	if (ap2 >= 0){
	    ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		    Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		    ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	    numExSphere++;
	    if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
	}
	else{
	    ErrWarnMessage(__LINE__, "Constraints::Delaunay:: no valid apex was found", 0);
	}

	//find the inclusion sphere (circumsphere of triangle n2-ap-n1)
	InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
	//this sphere could be too large or already contain other vertices
	//thus adding it in the containts set is useless
	if (InSphere[numInSphere].x[3] > 2){
	    //because the domain is scaled inside the unit box
	    n2 = ap;
	    continue;
	}

	aux_list[1] = ap;
	aux_list[2] = n1;
	aux_list[3] = n2;

	if (!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], nList, aux_list, Verts) ||
		!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[ap].connect, aux_list, Verts) ||
		!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n1].connect, aux_list, Verts) ||
		!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n2].connect, aux_list, Verts)){
	    //this circumsphere already contains other vertices and so it not a delaunay triangle
	    n2 = ap;

	    continue;
	}

	numInSphere++;
	if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
	n2 = ap;
    }
    return true;
}
bool Constraints::DelaunayNotConnected(int*nList, int*skipList, vert*Verts, bool loadSkipList)
{

    //nList is the sort list of neighbours around the a void but not a connected chain
    //i.e., last and first are not connected

    //1) find the exclusion spheres such that a new vertex should lie outside to maintain delaunayness of
    //untouched/surounding triangles
    //2) find the inclusion spheres such that a new vertex should lie within to create delaunay-accurate set
    //of triangles
    //skipList contians current vertices that should be negelected

    //start by loading skipList to aux_list
    //for insert, we always wanna consider the vertices in skipList (then set loadSkipList to false)
    if (loadSkipList){
	for (int i = 1; i <= skipList[0]; i++){
	    aux_list[i + 3] = skipList[i];
	}
	aux_list[0] = skipList[0] + 3;
    }
    else{
	aux_list[0] = 3;
    }

    if (nList[0] == 3){
	//this is easy because this result into one InSphere (that connect these three vertices)
	//and two ExSpheres (no looping needed)

	int ap(nList[1]), n1(nList[2]), n2(nList[3]), ap2;
	InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
	numInSphere++;


	//first exsphere
	ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;

	//second exsphere
	ap2 = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;
	return true;
    }





    for (int i = 1; i < nList[0]; i++){
	int ap = nList[i];
	int n1 = nList[i + 1];

	int ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);

	//find the exclusion sphere (circumsphere of triangle n2-ap-ap2)
	if (ap2 >= 0){
	    ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		    Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
		    ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	    numExSphere++;
	    if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
	}
	else{
	    return false;
	    ErrWarnMessage(__LINE__, "Constraints::DelaunayNotConnected:: no valid apex was found", 0);
	}

	if (nList[0] - i >=2){
	    int n2 = nList[i + 2];//two hubs away

	    //find the inclusion sphere (circumsphere of triangle n2-ap-n1)
	    InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
		    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
		    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
		    InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
	    //this sphere could be too large or already contain other vertices
	    //thus adding it in the containts set is useless
	    if (InSphere[numInSphere].x[3] > 2){
		//because the domain is scaled inside the unit box
		continue;
	    }

	    aux_list[1] = ap;
	    aux_list[2] = n1;
	    aux_list[3] = n2;

	    if (!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], nList, aux_list, Verts) ||
		    !IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[ap].connect, aux_list, Verts) ||
		    !IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n1].connect, aux_list, Verts) ||
		    !IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n2].connect, aux_list, Verts)){
		//this circumsphere already contains other vertices and so it not a delaunay triangle
		continue;
	    }
	    numInSphere++;
	    if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
	}
    }
    return true;
}
bool Constraints::DelaunayForcedNotConnected(int*nList, int removed, int*skipList, vert*Verts)
{
    //TODO add finding the circum spheres also (not well defined)
    //force nList to be not connected by removing one vertex from it
    //nList is the sort list of neighbours around the a void

    //1) find the exclusion spheres such that a new vertex should lie outside to maintain delaunayness of
    //untouched/surounding triangles


    //start by loading skipList to aux_list
    //for insert, we always wanna consider the vertices in skipList


    if (nList[0] == 3){
	return false;
	ErrWarnMessage(__LINE__, "Constraints::DelaunayForcedNotConnected:: have not considered this case yet", 0);
	//this is easy because this result into one InSphere (that connect these three vertices)
	//and two ExSpheres (no looping needed)
	/*int ap(nList[1]), n1(nList[2]), n2(nList[3]), ap2;
	  InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
	  Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
	  Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
	  InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
	  numInSphere++;


	//first exsphere
	ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
	Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
	Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
	ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;

	//second exsphere
	ap2 = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
	ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
	Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
	Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
	ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
	numExSphere++;
	skipList[0] -= 3;
	return;*/
    }


    int n1 = nList[nList[0]];
    for (int i = 1; i <= nList[0]; i++){
	int ap = nList[i];
	if (ap != removed && n1 != removed){
	    int ap2 = FindCommonElement_SkipList(Verts[ap].connect, Verts[n1].connect, skipList);

	    //find the exclusion sphere (circumsphere of triangle n2-ap-ap2)
	    if (ap2 >= 0){
		ExSphere[numExSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
			Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
			Verts[ap2].x[0], Verts[ap2].x[1], Verts[ap2].x[2],
			ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2]);
		numExSphere++;
		if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }

	    }
	    else{
		return false;
		ErrWarnMessage(__LINE__, "Constraints::Delaunay:: no valid apex was found", 0);
	    }
	}
	/*if (nList[0] - i > 2){
	  int n2 = nList[i + 2];//two hubs away

	//find the inclusion sphere (circumsphere of triangle n2-ap-n1)
	InSphere[numInSphere].x[3] = TriCircumcenter3d(Verts[ap].x[0], Verts[ap].x[1], Verts[ap].x[2],
	Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2],
	Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2],
	InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2]);
	//this sphere could be too large or already contain other vertices
	//thus adding it in the containts set is useless
	if (InSphere[numInSphere].x[3] > 2){
	//because the domain is scaled inside the unit box
	continue;
	}

	aux_list[1] = ap;
	aux_list[2] = n1;
	aux_list[3] = n2;

	if (!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], nList, aux_list, Verts) ||
	!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[ap].connect, aux_list, Verts) ||
	!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n1].connect, aux_list, Verts) ||
	!IsEmptySphere(InSphere[numInSphere].x[0], InSphere[numInSphere].x[1], InSphere[numInSphere].x[2], InSphere[numInSphere].x[3], Verts[n2].connect, aux_list, Verts)){
	//this circumsphere already contains other vertices and so it not a delaunay triangle
	continue;
	}
	numInSphere++;
	if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
	}*/
	n1 = ap;
    }
    return true;
}
bool Constraints::IsEmptySphere(double xc, double yc, double zc, double rc_2,
	int*list,
	int*skip_list,
	vert*Verts)
{
    for (int V = 1; V <= list[0]; V++){
	if (GetIndex(list[V], skip_list) >= 0){
	    continue;
	}
	double dist = Dist(Verts[list[V]].x[0], Verts[list[V]].x[1], Verts[list[V]].x[2], xc, yc, zc);
	if (dist<rc_2 - _tol*_tol){
	    return false;
	}
    }
    return true;
}
void Constraints::MinMaxAngle(int*nList, vert*Verts, double min_ang, double max_ang, double*void_vertex, bool nonObtuse, bool isConnected)
{
    //Find the inclusion and excluion region for all base angles such that a new vertex placed
    //in the void surounded by nList and connected to all vertices in nList will not violate
    //the max_ang and min_ang constraints (just the base angles)
    //in case of nonObtuse = true, we additionally add diameter spheres such that apex angle is less than 90.0

    int n2 = nList[nList[0]];

    for (int i = 1; i <= nList[0]; i++){
	int n1 = nList[i];

	if (i == 1 && !isConnected){
	    n2 = n1;
	    continue;
	}
	if (min_ang > _tol){
	    //*********** 1) Get min angle planes
	    //******* a) base angle
	    double rot_angle(90.0 - min_ang);
	    double x_n, y_n, z_n;

	    //normal to n1,n2, mid plane
	    Cross(Verts[n2].x[0] - Verts[n1].x[0], Verts[n2].x[1] - Verts[n1].x[1], Verts[n2].x[2] - Verts[n1].x[2],
		    void_vertex[0] - Verts[n1].x[0], void_vertex[1] - Verts[n1].x[1], void_vertex[2] - Verts[n1].x[2], x_n, y_n, z_n);
	    NormalizeVector(x_n, y_n, z_n);

	    //rotate n2 by rot_angle about axis normal to the plane containing n1,n2 and mid
	    //the axis passes through n1
	    double x1, y1, z1;
	    PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
		    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x1, y1, z1);
	    //To be improved/revised
	    //a naive way to check on the rotation direction
	    //by doing it again with negative sign
	    //and check distanced min_point
	    double x_sub, y_sub, z_sub;
	    rot_angle *= -1.0;
	    PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
		    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x_sub, y_sub, z_sub);

	    double len1(Dist(x_sub, y_sub, z_sub, void_vertex[0], void_vertex[1], void_vertex[2])), len2(Dist(x1, y1, z1, void_vertex[0], void_vertex[1], void_vertex[2]));
	    if (len1 > len2){
		x1 = x_sub;
		y1 = y_sub;
		z1 = z_sub;
	    }
	    else{
		rot_angle *= -1.0;
	    }


	    SinglePlane(x1 - Verts[n1].x[0],
		    y1 - Verts[n1].x[1],
		    z1 - Verts[n1].x[2],
		    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);


	    //same thing is done on the other end of the edge
	    Cross(Verts[n1].x[0] - Verts[n2].x[0], Verts[n1].x[1] - Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
		    void_vertex[0] - Verts[n2].x[0], void_vertex[1] - Verts[n2].x[1], void_vertex[2] - Verts[n2].x[2], x_n, y_n, z_n);
	    NormalizeVector(x_n, y_n, z_n);
	    double x2, y2, z2;
	    PerformRotation(Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], x_n, y_n, z_n,
		    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], rot_angle, x2, y2, z2);


	    SinglePlane(x2 - Verts[n2].x[0],
		    y2 - Verts[n2].x[1],
		    z2 - Verts[n2].x[2],
		    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2]);
	}


	if (max_ang > _tol){
	    //*********** 2) Get max angle planes
	    if (nonObtuse){
		//******* a) base angle
		if (!(*isObtuse)){
		    SinglePlane(Verts[n1].x[0] - Verts[n2].x[0],
			    Verts[n1].x[1] - Verts[n2].x[1],
			    Verts[n1].x[2] - Verts[n2].x[2],
			    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);

		    SinglePlane(Verts[n2].x[0] - Verts[n1].x[0],
			    Verts[n2].x[1] - Verts[n1].x[1],
			    Verts[n2].x[2] - Verts[n1].x[2],
			    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2]);
		}


		//******* b) apex angle
		ExSphere[numExSphere].x[0] = (Verts[n1].x[0] + Verts[n2].x[0]) / 2.0;
		ExSphere[numExSphere].x[1] = (Verts[n1].x[1] + Verts[n2].x[1]) / 2.0;
		ExSphere[numExSphere].x[2] = (Verts[n1].x[2] + Verts[n2].x[2]) / 2.0;
		ExSphere[numExSphere].x[3] = Dist(ExSphere[numExSphere].x[0], ExSphere[numExSphere].x[1], ExSphere[numExSphere].x[2], Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);

		numExSphere++;
		if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
	    }
	    else {
		// ??? always true?  why the isAcute check?
		if (true || !(*isAcute)){
		    //******* a) base angle
		    double rot_angle = abs(max_ang - 90.0);
		    double x_n, y_n, z_n;
		    //normal to n1,n2, mid plane
		    Cross(Verts[n2].x[0] - Verts[n1].x[0], Verts[n2].x[1] - Verts[n1].x[1], Verts[n2].x[2] - Verts[n1].x[2],
			    void_vertex[0] - Verts[n1].x[0], void_vertex[1] - Verts[n1].x[1], void_vertex[2] - Verts[n1].x[2], x_n, y_n, z_n);
		    NormalizeVector(x_n, y_n, z_n);

		    //rotate n2 by rot_angle about axis normal to the plane containing n1,n2 and mid
		    //the axis passes through n1
		    double x1, y1, z1;
		    PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
			    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x1, y1, z1);
		    //To be improved/revised
		    //a naive way to check on the rotation direction
		    //by doing it again with negative sign
		    //and check distanced min_point
		    double x_sub, y_sub, z_sub;
		    rot_angle *= -1.0;
		    PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n, y_n, z_n,
			    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], rot_angle, x_sub, y_sub, z_sub);

		    double len1 = Dist(x_sub, y_sub, z_sub, void_vertex[0], void_vertex[1], void_vertex[2]);
		    double len2 = Dist(x1, y1, z1, void_vertex[0], void_vertex[1], void_vertex[2]);
		    if (max_ang < 90.0 + _tol){ std::swap(len1, len2); }

		    if (len1 < len2){
			x1 = x_sub;
			y1 = y_sub;
			z1 = z_sub;
		    }
		    else {
			rot_angle *= -1.0;
		    }



		    SinglePlane(Verts[n1].x[0] - x1,
			    Verts[n1].x[1] - y1,
			    Verts[n1].x[2] - z1,
			    Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);


		    //same thing is done on the other end of the edge
		    Cross(Verts[n1].x[0] - Verts[n2].x[0], Verts[n1].x[1] - Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
			    void_vertex[0] - Verts[n2].x[0], void_vertex[1] - Verts[n2].x[1], void_vertex[2] - Verts[n2].x[2], x_n, y_n, z_n);
		    NormalizeVector(x_n, y_n, z_n);
		    double x2, y2, z2;
		    PerformRotation(Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], x_n, y_n, z_n,
			    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], rot_angle, x2, y2, z2);
		    SinglePlane(Verts[n2].x[0] - x2,
			    Verts[n2].x[1] - y2,
			    Verts[n2].x[2] - z2,
			    Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2]);
		}
	    }
	}

	n2 = n1;
    }

}
void Constraints::MinEdgeLength(int*nList, vert*Verts, double(*sizingfunc)(double xx, double yy, double zz))
{
    //min edge lenght constraints (to respect the sizing function)
    //it is based on the fact that new vertex in the void will be conntected to all
    //vertices in nList
    //sizingfunc should return the sizing function which is a protecteing sphere (ExSphere)
    //around each vertex in nList
    for (int i = 1; i <= nList[0]; i++){
	int n = nList[i];
	ExSphere[numExSphere].x[0] = Verts[n].x[0];
	ExSphere[numExSphere].x[1] = Verts[n].x[1];
	ExSphere[numExSphere].x[2] = Verts[n].x[2];
	ExSphere[numExSphere].x[3] = sizingfunc(Verts[n].x[0], Verts[n].x[1], Verts[n].x[2]);
	numExSphere++;
	if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
    }
}
void Constraints::MinEdgeLength(int*nList, vert*Verts, double r_min_2)
{
    //min edge lenght constraints (to respect the sizing function)
    //it is based on the fact that new vertex in the void will be conntected to all
    //vertices in nList
    //for unifrom sizing function, r_min_2 is square the sizing function
    for (int i = 1; i <= nList[0]; i++){
	int n = nList[i];
	ExSphere[numExSphere].x[0] = Verts[n].x[0];
	ExSphere[numExSphere].x[1] = Verts[n].x[1];
	ExSphere[numExSphere].x[2] = Verts[n].x[2];
	ExSphere[numExSphere].x[3] = r_min_2;

	numExSphere++;
	if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
    }
}
void Constraints::SinglePlane(double xn, double yn, double zn, double px, double py, double pz)
{
    //construct a single plance by specifiying the normal to the plane (xn, yn, zn)
    // and a point on it
    NormalizeVector(xn, yn, zn);

    Plane[numPlane].x[0] = xn;
    Plane[numPlane].x[1] = yn;
    Plane[numPlane].x[2] = zn;
    Plane[numPlane].x[3] = -1.0*(px * Plane[numPlane].x[0] + py * Plane[numPlane].x[1] + pz * Plane[numPlane].x[2]);

    numPlane++;
    if (numPlane == numPlane_size){ ExpandPlanes(numPlane_size, numPlane, Plane); }


}
void Constraints::Smoothness(int*nList, int*skipList, vert*Verts, double*void_vertex, double dev)
{
    //we are looking for two planes
    //one of them is define by the mid point between n1-n2
    //and has a normal of plane containing n1,n2 and ap but tilted an angle=
    //the other has the same point but the normal tilted with -ve of previous angle
    if (nList[0] == 3){
	//expand the skipList to include nList to get the correct apex
	skipList[++skipList[0]] = nList[1];
	skipList[++skipList[0]] = nList[2];
	skipList[++skipList[0]] = nList[3];
    }
    int n2 = nList[nList[0]];
    for (int i = 1; i <= nList[0]; i++){
	int n1 = nList[i];
	int ap = FindCommonElement_SkipList(Verts[n1].connect, Verts[n2].connect, skipList);
	double tri_tri_angle = TriTriNormalAngle(Verts[n1].x, Verts[n2].x, Verts[ap].x, void_vertex);
	double rot_ang = dev;

	//the plane normal vector
	double x_n, y_n, z_n;
	Cross(Verts[n1].x[0] - Verts[ap].x[0], Verts[n1].x[1] - Verts[ap].x[1], Verts[n1].x[2] - Verts[ap].x[2],
		Verts[n2].x[0] - Verts[ap].x[0], Verts[n2].x[1] - Verts[ap].x[1], Verts[n2].x[2] - Verts[ap].x[2],
		x_n, y_n, z_n);
	NormalizeVector(x_n, y_n, z_n);

	//point around which we will rotate the norm vector
	double x_mid = (Verts[n1].x[0] + Verts[n2].x[0]) / 2.0;
	double y_mid = (Verts[n1].x[1] + Verts[n2].x[1]) / 2.0;
	double z_mid = (Verts[n1].x[2] + Verts[n2].x[2]) / 2.0;

	//point on the vector
	double pv_x = x_n + x_mid;
	double pv_y = y_n + y_mid;
	double pv_z = z_n + z_mid;
	double pv_x_rot, pv_y_rot, pv_z_rot;

	//first plane
	PerformRotation(pv_x, pv_y, pv_z,
		Verts[n1].x[0] - Verts[n2].x[0], Verts[n1].x[1] - Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
		x_mid, y_mid, z_mid,
		rot_ang,
		pv_x_rot, pv_y_rot, pv_z_rot);

	Plane[numPlane].x[0] = pv_x_rot - x_mid;
	Plane[numPlane].x[1] = pv_y_rot - y_mid;
	Plane[numPlane].x[2] = pv_z_rot - z_mid;
	Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);

	if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
	    Plane[numPlane].x[0] *= -1.0;
	    Plane[numPlane].x[1] *= -1.0;
	    Plane[numPlane].x[2] *= -1.0;
	    Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);

	    if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
		ErrWarnMessage(__LINE__, "Constraints::Smoothness:: error(1)", 1);
	    }
	}

	numPlane++;
	if (numPlane == numPlane_size){ ExpandPlanes(numPlane_size, numPlane, Plane); }

	//second plane
	double rot_ang2 = -1.0*rot_ang;
	PerformRotation(pv_x, pv_y, pv_z,
		Verts[n1].x[0] - Verts[n2].x[0],Verts[n1].x[1] -Verts[n2].x[1], Verts[n1].x[2] - Verts[n2].x[2],
		x_mid, y_mid, z_mid,
		rot_ang2,
		pv_x_rot, pv_y_rot, pv_z_rot);

	Plane[numPlane].x[0] = pv_x_rot - x_mid;
	Plane[numPlane].x[1] = pv_y_rot - y_mid;
	Plane[numPlane].x[2] = pv_z_rot - z_mid;
	Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);


	if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
	    Plane[numPlane].x[0] *= -1.0;
	    Plane[numPlane].x[1] *= -1.0;
	    Plane[numPlane].x[2] *= -1.0;
	    Plane[numPlane].x[3] = -1.0*(x_mid*Plane[numPlane].x[0] + y_mid*Plane[numPlane].x[1] + z_mid*Plane[numPlane].x[2]);

	    if (Plane[numPlane].x[0] * void_vertex[0] + Plane[numPlane].x[1] * void_vertex[1] + Plane[numPlane].x[2] * void_vertex[2] + Plane[numPlane].x[3] > _tol){
		ErrWarnMessage(__LINE__, "Constraints::Smoothness:: error(2)", 1);
	    }
	}

	numPlane++;
	if (numPlane == numPlane_size){ ExpandPlanes(numPlane_size, numPlane, Plane); }

	n2 = n1;
    }
    if (nList[0] == 3){
	skipList[0] -= 3;
    }
}
void Constraints::Maximality()
{
    //TODO
}
void Constraints::AttractorInSphere(vert*Verts, int ip, double*void_vertex)
{
    //construct an inclusion sphere around void_vertex such that wehn ip is relocated
    //it is guarantteed to be closer to the void_vertex
    InSphere[numInSphere].x[0] = void_vertex[0];
    InSphere[numInSphere].x[1] = void_vertex[1];
    InSphere[numInSphere].x[2] = void_vertex[2];
    InSphere[numInSphere].x[3] = Dist(Verts[ip].x[0], Verts[ip].x[1], Verts[ip].x[2], void_vertex[0], void_vertex[1], void_vertex[2]);
    numInSphere++;
    if (numInSphere == numInSphere_size){ ExpandSpheres(numInSphere_size, numInSphere, InSphere); }
}
void Constraints::RepellerExSphere(vert*Verts, int ip, double*void_vertex)
{
    //construct an exclusion sphere around void_vertex such that wehn ip is relocated
    //it is guarantteed to be further away from the void_vertex
    ExSphere[numExSphere].x[0] = void_vertex[0];
    ExSphere[numExSphere].x[1] = void_vertex[1];
    ExSphere[numExSphere].x[2] = void_vertex[2];
    ExSphere[numExSphere].x[3] = Dist(Verts[ip].x[0], Verts[ip].x[1], Verts[ip].x[2], void_vertex[0], void_vertex[1], void_vertex[2]);
    numExSphere++;
    if (numExSphere == numExSphere_size){ ExpandSpheres(numExSphere_size, numExSphere, ExSphere); }
}
void Constraints::Direction(int*ip, int*nList, vert*Verts)
{
    //this builds a bounding polygon to prevent a newly relocated or created
    //vertex from being created in a region that might create a tangled mesh
    //while preserving all other quality matric

    //loop around nList
    //for each edge in nList starting by n1-n2
    //find the common vertex of n1 and n2 in ip --> shrd
    //rotate the line between mid and n1 by 90 in direction closer to shrd
    //construct the plane with normal as the rotate line
    //and the point mid(n1-n2) as a point in it

    //nList could be connected or not connected (does not make difference)

    int n2 = nList[nList[0]];
    for (int i = 1; i <= nList[0]; i++){
	int n1 = nList[i];

	//find the share vertex
	if (FindCommonElements(Verts[n1].connect, Verts[n2].connect, aux_list)){
	    int id = -1;
	    for (int j = 1; j <= aux_list[0]; j++){
		int sh = aux_list[j];
		id = (ip != NULL) ? GetIndex(aux_list[j], ip) : GetIndex(aux_list[j], nList); //for injection, the shared vertex is in nList

		/*if (j == aux_list[0] && id < 0 && nList[0] == 3 && ip[0] == 1){
		//special case of nList not connected with three vertices only
		id = 0;
		sh = ip[1];
		}*/
		if (id >= 0){

		    double xmid(0.5*(Verts[n1].x[0] + Verts[n2].x[0])),
			   ymid(0.5*(Verts[n1].x[1] + Verts[n2].x[1])),
			   zmid(0.5*(Verts[n1].x[2] + Verts[n2].x[2]));

		    //normal to n1,n2, sh
		    double x_n12, y_n12, z_n12;
		    Cross(Verts[n2].x[0] - Verts[n1].x[0], Verts[n2].x[1] - Verts[n1].x[1], Verts[n2].x[2] - Verts[n1].x[2],
			    Verts[sh].x[0] - Verts[n1].x[0], Verts[sh].x[1] - Verts[n1].x[1], Verts[sh].x[2] - Verts[n1].x[2], x_n12, y_n12, z_n12);
		    NormalizeVector(x_n12, y_n12, z_n12);

		    double x1, y1, z1;
		    //PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
		    //	            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], 90.0, x1, y1, z1);
		    PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
			    xmid, ymid, zmid, 90.0, x1, y1, z1);
		    double x_sub, y_sub, z_sub;
		    //PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
		    //	            Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2], -90.0, x_sub, y_sub, z_sub);
		    PerformRotation(Verts[n2].x[0], Verts[n2].x[1], Verts[n2].x[2], x_n12, y_n12, z_n12,
			    xmid, ymid, zmid, -90.0, x_sub, y_sub, z_sub);

		    double len1(Dist(x_sub, y_sub, z_sub, Verts[sh].x[0], Verts[sh].x[1], Verts[sh].x[2])), len2(Dist(x1, y1, z1, Verts[sh].x[0], Verts[sh].x[1], Verts[sh].x[2]));
		    if (len1 > len2){
			x1 = x_sub;
			y1 = y_sub;
			z1 = z_sub;
		    }

		    //normal
		    //double xn(x1 - Verts[n1].x[0]),
		    //	     yn(y1 - Verts[n1].x[1]),
		    //	     zn(z1 - Verts[n1].x[2]);

		    double xn(x1 - xmid),
			   yn(y1 - ymid),
			   zn(z1 - zmid);
		    SinglePlane(xn, yn, zn, xmid, ymid, zmid);
		    //SinglePlane(xn, yn, zn, Verts[n1].x[0], Verts[n1].x[1], Verts[n1].x[2]);
		    break;
		}
	    }
	    if (id < 0 && nList[0] != 3){
		ErrWarnMessage(__LINE__, "Constraints::Direction edge has incorrect shared vertex!!!", 0);
	    }
	}
	else{
	    ErrWarnMessage(__LINE__, "Constraints::Direction edge does not have shared vertex!!!", 0);
	}
	n2 = n1;
    }
}
//*** Checking
bool Constraints::OverlappingInSpheres()
{
    //find if all InSpheres overlaps
    //if two don't overlap, return false
    for (int i = 0; i < numInSphere - 1; i++){
	for (int j = i + 1; j < numInSphere; j++){
	    if (Dist(InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2], InSphere[j].x[0], InSphere[j].x[1], InSphere[j].x[2]) *_tol_sq_circ >
		    InSphere[i].x[3] + InSphere[j].x[3] + 2.0*sqrt(InSphere[i].x[3] * InSphere[j].x[3])){
		return false;
	    }
	}
    }

    return true;
}
bool Constraints::InsideFeasibleRegion_Vertex(double xx, double yy, double zz)
{
    //check if (xx,yy,zz) is insdie all the feasible region
    //i.e., inside all inclusion regions and outside all exclusion regions
    if (!IsInsideBoundingBox(xx, yy, zz)){
	return false;
    }

    for (int i = 0; i < numExSphere; i++){
	double dist = Dist(xx, yy, zz, ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);
	if (dist < ExSphere[i].x[3] + _tol_sq){ return false; }
    }

    for (int i = 0; i < numInSphere; i++){
	double dist = Dist(xx, yy, zz, InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);
	if (dist > InSphere[i].x[3] - _tol_sq){ return false; }
    }

    for (int i = 0; i < numPlane; i++){
	if (xx*Plane[i].x[0] + yy*Plane[i].x[1] + zz*Plane[i].x[2] + Plane[i].x[3]>-_tol){
	    return false;
	}
    }

    return true;
}
bool Constraints::InsideFeasibleRegion_Triangle(double*tri)
{
    //tri is the coordinates of the triangle
    //0-2 is first vertex
    //3-5 is second vertex
    //6-8 is third vertex

    for (int i = 0; i < numExSphere; i++){
	double dist1 = Dist(tri[0], tri[1], tri[2], ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);
	double dist2 = Dist(tri[3], tri[4], tri[5], ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);
	double dist3 = Dist(tri[6], tri[7], tri[8], ExSphere[i].x[0], ExSphere[i].x[1], ExSphere[i].x[2]);

	if (dist1<ExSphere[i].x[3] * _tol_sq_circ && dist2<ExSphere[i].x[3] * _tol_sq_circ && dist3<ExSphere[i].x[3] * _tol_sq_circ){
	    //the whole triangle is in one of the exclusion sphere
	    return false;
	}
    }

    for (int i = 0; i < numInSphere; i++){
	double dist1 = Dist(tri[0], tri[1], tri[2], InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);
	double dist2 = Dist(tri[3], tri[4], tri[5], InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);
	double dist3 = Dist(tri[6], tri[7], tri[8], InSphere[i].x[0], InSphere[i].x[1], InSphere[i].x[2]);

	if (dist1 > InSphere[i].x[3] * _tol_sq_circ && dist2 > InSphere[i].x[3] * _tol_sq_circ && dist3 > InSphere[i].x[3] * _tol_sq_circ){
	    //the whole triangle is in one of the exclusion sphere
	    return false;
	}
    }


    for (int i = 0; i < numPlane; i++){
	if (tri[0] * Plane[i].x[0] + tri[1] * Plane[i].x[1] + tri[2] * Plane[i].x[2] + Plane[i].x[3] >-_tol &&
		tri[3] * Plane[i].x[0] + tri[4] * Plane[i].x[1] + tri[5] * Plane[i].x[2] + Plane[i].x[3] >-_tol &&
		tri[6] * Plane[i].x[0] + tri[7] * Plane[i].x[1] + tri[8] * Plane[i].x[2] + Plane[i].x[3] >-_tol ){

	    return false;
	}
    }
    return true;

}
bool Constraints::IsInsideBoundingBox(double xx, double yy, double zz)
{
    //check if (xx,yy,zz) is inside the bounding box myBox
    if (abs(myBox.xmin - DBL_MAX) < _tol){
	//if the box did not set, we assume it is okay
	return true;
    }
    else{
	if (xx < myBox.xmin || xx > myBox.xmax ||
		yy < myBox.ymin || yy > myBox.ymax ||
		zz < myBox.zmin || xx > myBox.zmax ){
	    return false;
	}
    }
    return true;
}


void Constraints::ExpandSpheres(int&currentSize,int currentNumSphere, sphere*&mySpheres)
{
    currentSize *= 2;
    sphere*newSpheres = new sphere[currentSize];
    for (int i = 0; i < currentNumSphere; i++){
	for (int j = 0; j < 4; j++){
	    newSpheres[i].x[j] = mySpheres[i].x[j];
	}
    }

    delete[]mySpheres;

    mySpheres = new sphere[currentSize];
    for (int i = 0; i < currentNumSphere; i++){
	for (int j = 0; j < 4; j++){
	    newSpheres[i].x[j] = mySpheres[i].x[j];
	}
    }
    currentSize--;
}
void Constraints::ExpandPlanes(int&currentSize, int currentNumPlane, plane*&myPlanes)
{
    currentSize *= 2;
    sphere*newPlanes = new sphere[currentSize];
    for (int i = 0; i < currentNumPlane; i++){
	for (int j = 0; j < 4; j++){
	    newPlanes[i].x[j] = myPlanes[i].x[j];
	}
    }

    delete[]myPlanes;

    myPlanes = new plane[currentSize];
    for (int i = 0; i < currentNumPlane; i++){
	for (int j = 0; j < 4; j++){
	    newPlanes[i].x[j] = myPlanes[i].x[j];
	}
    }
    currentSize--;
}

//*** Debug
void Constraints::SpheresDrawing(std::string filename, int num, sphere *Sphere)
{
    double **drawpsheres = new double*[num];
    for (int i = 0; i < num; i++){
	drawpsheres[i] = new double[4];
	drawpsheres[i][0] = Sphere[i].x[0];
	drawpsheres[i][1] = Sphere[i].x[1];
	drawpsheres[i][2] = Sphere[i].x[2];
	drawpsheres[i][3] = Sphere[i].x[3];
    }

    //DrawManySpheres(filename, num, drawpsheres, 1);
    for (int i = 0; i < num; i++){
	delete[] drawpsheres[i];
    }
    free(drawpsheres);
}
void Constraints::DrawInSpheres()
{
    ///SpheresDrawing("debug_out/InSphere.obj",numInSphere, InSphere);
}
void Constraints::DrawExSpheres()
{
    //SpheresDrawing("debug_out/ExSphere.obj", numExSphere, ExSphere);
}

double OpimizerFunc_CenterAngle(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);
double OpimizerFunc_SideAngle(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz);
double OpimizerFunc_Closer(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);
double OpimizerFunc_Further(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz);

Operator::Operator(vert*Vert_org, bool *isO, bool *isA):
    Vert_org(Vert_org)
{
    myConstraints.isObtuse = isO;
    myConstraints.isAcute = isA;
    active_pool = new int *[MAXPOOL];
    tri_pool = new int *[MAXPOOL];
    tmp_active_pool = new int*[MAXPOOL];
    for (int i = 0; i < MAXPOOL; i++){
	active_pool[i] = new int[3];
	tri_pool[i] = new int[3];
	tmp_active_pool[i] = new int[3];
    }

    _tar = new double*[4];
    _tar[0] = new double[9];
    _tar[1] = new double[9];
    _tar[2] = new double[9];
    _tar[3] = new double[9];


    srand(time(NULL));
    myRandNum.InitiateRndNum(rand());

}

Operator::~Operator()
{

}

///**************Special
void Operator::TriValentRemoval(int&numVert, vert*Verts)
{
    //remove all trivalent regardless to their effect on the quality
    //do it multiple times
    //everything can be fixed later
    //better to call this before applying the operators
    //some operators (injection) quit when it finds tri-valents on a patch since it will mess up the toplogical correctness
    //other operators can tolerate tri-valents vertices as long as topological correctness holds
    while (true){
	bool again(false);
	for (int i = numVert; i >= 0; i--){
	    if (Verts[i].connect[0] == 3){
		RemoveVertex(i, numVert, Verts);
		again = true;
	    }
	}
	if (!again){ break; }
    }
}



///**************Operators
bool Operator::Relocation(int ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer)
{
    //TODO pass the optimizer function also
    //ip = the vertex to relocate
    //numVert =  total num of vertices in mesh
    //Vert =  the mesh to update



    myConstraints.Reset(Verts, NULL);

    skipList[0] = 1;
    skipList[1] = ip;

    if (constraints.isDelaunay){
	//Delaunay constraints
	myConstraints.Delaunay(Verts[ip].connect, skipList, Verts, true,1);
	if (!myConstraints.OverlappingInSpheres()){
	    //empty inclusion region
	    return false;
	}
    }

    myConstraints.Direction(skipList, Verts[ip].connect, Verts);

    if (constraints.isMinAngle || constraints.isMaxAngle){
	//Min max angle
	myConstraints.MinMaxAngle(Verts[ip].connect, Verts, constraints.MinAngle, constraints.MaxAngle, Verts[ip].x, constraints.isNonobtuse,1);
    }

    if (constraints.isEdgeLen){
	//Min edge length
	myConstraints.MinEdgeLength(Verts[ip].connect, Verts, constraints.MinEdgeLength_sq);
    }

    if (constraints.isSmooth){
	//Smoothness
	myConstraints.Smoothness(Verts[ip].connect, skipList, Verts, Verts[ip].x, constraints.dev);
    }


    StartActivePool(closedtSurfaceID, numSurfaceLayer);
    if (Sampler(Verts, Verts[ip].connect, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
	Mover(ip, Verts);
	return true;
    }
    return false;

}
bool Operator::Ejection(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool att)
{
    //TODO pass the optimizer function also
    //ip = pointer to vertices to eject
    //numVert =  total num of vertices in mesh
    //Vert =  the mesh to update

#ifdef DEBUGGING
    if (ip[0] != 2 && ip[0] != 3){
	ErrWarnMessage(__LINE__, "Operator::Ejection:: can only ejection two or three vertices", 0);
    }
#endif

    //sanity check
    //make sure vertices in ip are not connected to the same tri valent node

    skipList[0] = 0;//use this list to store the tri valent nodes connected to any of the to-be ejected vertices
    TriValent(ip[1], Verts, skipList, ip[2], (ip[0] == 3) ? ip[3] : INT_MAX);
    TriValent(ip[2], Verts, skipList, ip[1], (ip[0] == 3) ? ip[3] : INT_MAX);
    if (ip[0] == 3){
	TriValent(ip[3], Verts, skipList, ip[1], ip[2]);
    }
    if (FindDuplication(skipList)){ return false; }

    //get mid vertex and copy skipList
    skipList[0] = ip[0];
    void_vertex[0] = void_vertex[1] = void_vertex[2] = 0;
    for (int i = 1; i <= ip[0]; i++){
	skipList[i] = ip[i];
	void_vertex[0] += Verts[ip[i]].x[0];
	void_vertex[1] += Verts[ip[i]].x[1];
	void_vertex[2] += Verts[ip[i]].x[2];
    }
    void_vertex[0] /= double(ip[0]); void_vertex[1] /= double(ip[0]); void_vertex[2] /= double(ip[0]);


    //get the sort neighbout list (unduplicated list of vertices connected to *ip and sorted)
    mynList[0] = 0;
    if (ip[0] == 2){
	if (Verts[ip[1]].connect[0] == 3){
	    GetListSkipVertex(ip[2], Verts, ip[1],INT_MAX, mynList);
	}
	else if (Verts[ip[2]].connect[0] == 3){
	    GetListSkipVertex(ip[1], Verts, ip[2], INT_MAX, mynList);
	}
	else{
	    GetEdgeSortedNeighbourList(ip[1], ip[2], Verts, mynList);
	    if (mynList[0] < 0 || FindDuplication(mynList)){
		return false;
	    }
	}
    }
    else{
	//TODO get the list correct when ejecting three vertices
	//take care of special cases of tri-valent nodes

	//it is not possible to have two tri-valent node conneced in a watertight manifold
	if (Verts[ip[1]].connect[0] == 3 || Verts[ip[2]].connect[0] == 3 || Verts[ip[3]].connect[0] == 3){
	    return false;
	}
	else {
	    GetFaceSortedNeighbourList(ip[1], ip[2], ip[3], Verts, mynList);
	    if (mynList[0] < 0 || FindDuplication(mynList)){
		return false;
	    }
	}
    }
    if (mynList[0] == 0){
	//could not get the neighbours right
	return false;
    }

    ///******* Attractor
    if (att){
	SetAttractorRepeller(Verts);
	AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, true);
    }

    myConstraints.Reset(Verts, NULL);
    myConstraints.Direction(ip, mynList, Verts);

    if (constraints.isDelaunay){
	//Delaunay constraints
	myConstraints.Delaunay(mynList, skipList, Verts, true, 1);
	if (!myConstraints.OverlappingInSpheres()){
	    //empty inclusion region
	    if (att){
		RevertAttractorRepeller(Verts);
	    }
	    return false;
	}
    }

    if (constraints.isMinAngle || constraints.isMaxAngle){
	//Min max angle
	myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse,1);
    }

    if (constraints.isEdgeLen){
	//Min edge length
	myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
    }
    if (constraints.isSmooth){
	//Smoothness
	myConstraints.Smoothness(mynList, ip, Verts, void_vertex, constraints.dev);
    }

    StartActivePool(closedtSurfaceID, numSurfaceLayer);
    if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
	//for update, call EdgeCollapse if it is two vertices or
	//FaceCollapse if it is three vertices
	if (ip[0] == 2){
	    EdgeCollapse(ip[1],ip[2],numVert,Verts);
	}
	else {
	    FaceCollapse(ip[1], ip[2], ip[3], numVert, Verts);
	}
	return true;
    }

    if (att){
	RevertAttractorRepeller(Verts);
    }

    return false;
}
bool Operator::AggressiveInjection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj)
{
    //TODO pass the optimizer function also
    //TODO check that three vertices in ip forms a triangle

    //ip = pointer to vertices to eject (should be three vertices)
    //numVert =  total num of vertices in mesh
    //Vert =  the mesh to update

#ifdef DEBUGGING
    if (ip[0] != 3){
	ErrWarnMessage(__LINE__, "Operator::Injection:: Invalud input. Correct input is three vertices (triangle)", 0);
    }
#endif

    if (Verts[ip[1]].connect[0] == 3 || Verts[ip[2]].connect[0] == 3 || Verts[ip[3]].connect[0] == 3){
	//nop, we don't do this
	return false;
    }


    skipList[0] = ip[0];
    void_vertex[0] = void_vertex[1] = void_vertex[2] = 0.0;
    for (int i = 1; i <= ip[0]; i++){
	skipList[i] = ip[i];
	void_vertex[0] += Verts[ip[i]].x[0];
	void_vertex[1] += Verts[ip[i]].x[1];
	void_vertex[2] += Verts[ip[i]].x[2];
    }
    void_vertex[0] /= double(ip[0]);
    void_vertex[1] /= double(ip[0]);
    void_vertex[2] /= double(ip[0]);

    //find the correct shared vertex between each two consecutive vertices in ip
    apex[0] = ip[0];
    int i_skip = -1;
    for (int i = 1; i <= ip[0]; i++){
	int iq = (i == ip[0]) ? ip[1] : ip[i + 1];
	//apex[i] = FindCommonElement_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip);

	if (FindCommonElements_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip, temp_arr)){
	    if (temp_arr[0] == 1){
		apex[i] = temp_arr[1];
	    }
	    else if (temp_arr[0] == 2){

		if (Verts[temp_arr[1]].connect[0] == 3 || Verts[temp_arr[2]].connect[0] == 3){
		    //get the trivalen one
		    apex[i] = (Verts[temp_arr[1]].connect[0] == 3) ? temp_arr[1] : temp_arr[2];
		}
		else{

		    //in this case we seek the vertex that is not connected to the other apex
		    //since other apex's are not discovered yet, we skip this one and do it after getting
		    //other apexs. There should not be a way that there is more than one apex that is such problomatic
		    //but we check on this anyways
#ifdef DEBUGGING
		    if (i_skip >= 0){
			ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
		    }
#endif
		    i_skip = i;
		}
	    }
	    else{
#ifdef DEBUGGING
		ErrWarnMessage(__LINE__, "Operator::Injection:: 2 can not get apex", 0);
#endif
		return false;
	    }

	}
	else{
	    return false;
#ifdef DEBUGGING
	    ErrWarnMessage(__LINE__, "Operator::Injection:: 1 can not get apex", 0);
#endif
	    return false;
	}

    }

    if (i_skip >= 0){
	//we have on problomatic apex
	//at least this apex should not be connected to
	int iq = (i_skip == ip[0]) ? ip[1] : ip[i_skip + 1];

	FindCommonElements_SkipList(Verts[ip[i_skip]].connect, Verts[iq].connect, ip, temp_arr);

	bool one_shared(false), two_shared(false);
	for (int j = 1; j <= 3; j++){
	    if (j == i_skip){ continue; }
	    if (GetIndex(temp_arr[1], Verts[apex[j]].connect) >= 0){ one_shared = true; }
	    if (GetIndex(temp_arr[2], Verts[apex[j]].connect) >= 0){ two_shared = true; }
	}

#ifdef DEBUGGING
	if ((one_shared&&two_shared) || (!one_shared&&!two_shared)){
	    //means both of the temp_arr[1] and temp_arr[2] (candidate apex)
	    //are either shared with other apex's and not shared at all
	    //have not considered this case yet
	    ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
	}
#endif

	apex[i_skip] = (two_shared) ? temp_arr[1] : temp_arr[2];
    }

    if (Verts[apex[1]].connect[0] == 3 || Verts[apex[2]].connect[0] == 3 || Verts[apex[3]].connect[0] == 3){
	//nop, we don't do this either
	return false;
    }




    //populate the void corner list
    mynList[0] = 6;
    mynList[1] = ip[1];
    mynList[2] = apex[1];
    mynList[3] = ip[2];
    mynList[4] = apex[2];
    mynList[5] = ip[3];
    mynList[6] = apex[3];


    //*******1) Destory the whole triangle as if we refining it
    //and removing the edges ip[1]-ip[2], ip[2]-ip[3] and ip[3]-ip[1]
    //thus creating a void surrounded by six vertices
    //only if non of the three vertices in ip are 4-valent
    //because when we insert a node in the middle in this fashion,
    //we effectively, reduce the valence of each node in ip by one

    if (Verts[ip[1]].connect[0] > 4 && Verts[ip[2]].connect[0] > 4 && Verts[ip[3]].connect[0] > 4 /*&& InspectFeasibleRegion(mynList, Verts)*/){

	if (inj){
	    SetAttractorRepeller(Verts);
	    AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
	}

	myConstraints.Reset(Verts, NULL);

	myConstraints.Direction(NULL, mynList, Verts);

	if (constraints.isDelaunay){
	    //Delaunay constraints
	    myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
	}

	if (myConstraints.OverlappingInSpheres()){
	    if (constraints.isMinAngle || constraints.isMaxAngle){
		//Min max angle
		myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
	    }

	    if (constraints.isEdgeLen){
		//Min edge length
		myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
	    }
	    if (constraints.isSmooth){
		//Smoothness
		myConstraints.Smoothness(mynList, ip, Verts, void_vertex, constraints.dev);
	    }

	    StartActivePool(closedtSurfaceID, numSurfaceLayer);
	    if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
		//load the six edges to break
		//ip[1]-ip[2]
		//ip[1]-ip[3]
		//ip[2]-ip[3]
		EdgeToBreak[0] = 3;
		EdgeToBreak[1] = ip[1];
		EdgeToBreak[2] = ip[2];
		EdgeToBreak[3] = ip[1];
		EdgeToBreak[4] = ip[3];
		EdgeToBreak[5] = ip[2];
		EdgeToBreak[6] = ip[3];
		Inserter(numVert, Verts);
		return true;
	    }
	}
	if (inj){
	    RevertAttractorRepeller(Verts);
	}
    }


    //*******2) Destory one edge of the triangle and another one
    //there is 9 possibilities in near, we try them all :D

    int ip3 = ip[ip[0]];
    for (int i = 1; i <= ip[0]; i++){
	int j = (i == ip[0]) ? 1 : i + 1;

	int ip1 = ip[i];
	int ip2 = ip[j];
	int iapex = apex[i];
	int japex = apex[j];

	//we take the void_vertex as the avergae of two points
	//point one is the mid-point between ip1-ip2
	//point one is the mid-point between ip2-ip3
	void_vertex[0] = (Verts[ip1].x[0] + 2.0*Verts[ip2].x[0] + Verts[ip3].x[0]) / 4.0;
	void_vertex[1] = (Verts[ip1].x[1] + 2.0*Verts[ip2].x[1] + Verts[ip3].x[1]) / 4.0;
	void_vertex[2] = (Verts[ip1].x[2] + 2.0*Verts[ip2].x[2] + Verts[ip3].x[2]) / 4.0;


	//#1 first try to destroy the edges ip1-ip2 ip2-ip3
	if (Verts[ip2].connect[0] > 4){ //this prevents ip2 from being tri-valent
					//(because we remove two verices connected to it and add one, effectively reduce valence by 1)

	    mynList[0] = 5;
	    mynList[1] = ip1;
	    mynList[2] = iapex;
	    mynList[3] = ip2;
	    mynList[4] = japex;
	    mynList[5] = ip3;

	    skipList[0] = 3;

	    if (inj){
		SetAttractorRepeller(Verts);
		AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
	    }

	    myConstraints.Reset(Verts, NULL);
	    myConstraints.Direction(NULL, mynList, Verts);

	    if (constraints.isDelaunay){
		//Delaunay constraints
		myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
	    }

	    if (myConstraints.OverlappingInSpheres()){
		if (constraints.isMinAngle || constraints.isMaxAngle){
		    //Min max angle
		    myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
		}

		if (constraints.isEdgeLen){
		    //Min edge length
		    myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
		}

		if (constraints.isSmooth){
		    //Smoothness
		    myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
		}


		StartActivePool(closedtSurfaceID, numSurfaceLayer);
		if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
		    //load edge to break
		    //ip1-ip2
		    //ip2-ip3
		    EdgeToBreak[0] = 2;
		    EdgeToBreak[1] = ip1;
		    EdgeToBreak[2] = ip2;
		    EdgeToBreak[3] = ip2;
		    EdgeToBreak[4] = ip3;
		    Inserter(numVert, Verts);
		    return true;
		}
	    }

	    if (inj){
		RevertAttractorRepeller(Verts);
	    }
	}



	//#2 second try to destroy the edges ip1-ip2 with another edge connected to on of the apex

	for (int id = 1; id <= 2; id++){
	    if (id == 2){
		std::swap(ip1, ip2);//trust me you need this.
	    }

	    if (Verts[ip1].connect[0] > 4){//this prevents ip1 from being tri-valent
					   //(because we remove to verices connected to it and add one, effectively reduce valence by 1)

					   //the two edges to destroy are ip1-ip2  and iapex1-ip1
		int extra_apex = FindCommonElement_SkipList(Verts[ip1].connect, Verts[iapex].connect, ip);

		if (Verts[extra_apex].connect[0] > 3){//if it is 3-valent, then the void is ambigious
						      //the only solution would be creating 4-valent vertex which is equally bad for most cases
		    mynList[0] = 5;
		    mynList[1] = ip1;
		    mynList[2] = extra_apex;
		    mynList[3] = iapex;
		    mynList[4] = ip2;
		    mynList[5] = ip3;

		    //here we take the void_vertex as the avergae of two points
		    //point one is the mid-point between ip1-iapex
		    //point one is the mid-point between ip1-ip2
		    void_vertex[0] = (2.0*Verts[ip1].x[0] + Verts[iapex].x[0] + Verts[ip2].x[0]) / 4.0;
		    void_vertex[1] = (2.0*Verts[ip1].x[1] + Verts[iapex].x[1] + Verts[ip2].x[1]) / 4.0;
		    void_vertex[2] = (2.0*Verts[ip1].x[2] + Verts[iapex].x[2] + Verts[ip2].x[2]) / 4.0;

		    skipList[0] = 4;
		    skipList[skipList[0]] = iapex;

		    if (inj){
			SetAttractorRepeller(Verts);
			AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
		    }

		    myConstraints.Reset(Verts, NULL);
		    myConstraints.Direction(NULL, mynList, Verts);

		    if (constraints.isDelaunay){
			//Delaunay constraints
			myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
		    }

		    if (myConstraints.OverlappingInSpheres()){
			if (constraints.isMinAngle || constraints.isMaxAngle){
			    //Min max angle
			    myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
			}

			if (constraints.isEdgeLen){
			    //Min edge length
			    myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
			}

			if (constraints.isSmooth){
			    //Smoothness
			    myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
			}


			StartActivePool(closedtSurfaceID, numSurfaceLayer);
			if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
			    EdgeToBreak[0] = 2;
			    EdgeToBreak[1] = ip1;
			    EdgeToBreak[2] = ip2;
			    EdgeToBreak[3] = ip1;
			    EdgeToBreak[4] = iapex;
			    Inserter(numVert, Verts);
			    return true;
			}
		    }

		    if (inj){
			RevertAttractorRepeller(Verts);
		    }
		}
	    }
	}

	ip3 = ip2;//not ip1 because we swap it

    }



    return false;
}
bool Operator::Injection(int*ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool inj)
{
    //TODO pass the optimizer function also
    //TODO check that three vertices in ip forms a triangle

    //ip = pointer to vertices to eject (should be three vertices)
    //numVert =  total num of vertices in mesh
    //Vert =  the mesh to update

#ifdef DEBUGGING
    if (ip[0] != 3){
	ErrWarnMessage(__LINE__, "Operator::Injection:: Invalud input. Correct input is three vertices (triangle)", 0);
    }
#endif

    if (Verts[ip[1]].connect[0] == 3 || Verts[ip[2]].connect[0] == 3 || Verts[ip[3]].connect[0] == 3){
	//nop, we don't do this
	return false;
    }


    skipList[0] = ip[0];
    void_vertex[0] = void_vertex[1] = void_vertex[2] = 0.0;
    for (int i = 1; i <= ip[0]; i++){
	skipList[i] = ip[i];
	void_vertex[0] += Verts[ip[i]].x[0];
	void_vertex[1] += Verts[ip[i]].x[1];
	void_vertex[2] += Verts[ip[i]].x[2];
    }
    void_vertex[0] /= double(ip[0]);
    void_vertex[1] /= double(ip[0]);
    void_vertex[2] /= double(ip[0]);

    //find the correct shared vertex between each two consecutive vertices in ip
    apex[0] = ip[0];
    int i_skip = -1;
    for (int i = 1; i <= ip[0]; i++){
	int iq = (i == ip[0]) ? ip[1] : ip[i + 1];
	//apex[i] = FindCommonElement_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip);

	if (FindCommonElements_SkipList(Verts[ip[i]].connect, Verts[iq].connect, ip, temp_arr)){
	    if (temp_arr[0] == 1){
		apex[i] = temp_arr[1];
	    }
	    else if (temp_arr[0] == 2){

		if (Verts[temp_arr[1]].connect[0] == 3 || Verts[temp_arr[2]].connect[0] == 3){
		    //get the trivalen one
		    apex[i] = (Verts[temp_arr[1]].connect[0] == 3) ? temp_arr[1] : temp_arr[2];
		}
		else{

		    //in this case we seek the vertex that is not connected to the other apex
		    //since other apex's are not discovered yet, we skip this one and do it after getting
		    //other apexs. There should not be a way that there is more than one apex that is such problomatic
		    //but we check on this anyways
#ifdef DEBUGGING
		    if (i_skip >= 0){
			ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
		    }
#endif
		    i_skip = i;
		}
	    }
	    else{
		return false;
		ErrWarnMessage(__LINE__, "Operator::Injection:: 2 can not get apex", 0);
	    }

	}
	else{
	    return false;
	    ErrWarnMessage(__LINE__, "Operator::Injection:: 1 can not get apex", 0);
	}

    }

    if (i_skip >= 0){
	//we have on problomatic apex
	//at least this apex should not be connected to
	int iq = (i_skip == ip[0]) ? ip[1] : ip[i_skip + 1];

	FindCommonElements_SkipList(Verts[ip[i_skip]].connect, Verts[iq].connect, ip, temp_arr);

	bool one_shared(false), two_shared(false);
	for (int j = 1; j <= 3; j++){
	    if (j == i_skip){ continue; }
	    if (GetIndex(temp_arr[1], Verts[apex[j]].connect) >= 0){ one_shared = true; }
	    if (GetIndex(temp_arr[2], Verts[apex[j]].connect) >= 0){ two_shared = true; }
	}

#ifdef DEBUGGING
	if ((one_shared&&two_shared) || (!one_shared&&!two_shared)){
	    //means both of the temp_arr[1] and temp_arr[2] (candidate apex)
	    //are either shared with other apex's and not shared at all
	    //have not considered this case yet
	    ErrWarnMessage(__LINE__, "Operator::Injection:: not considered. See comment", 0);
	}
#endif

	apex[i_skip] = (two_shared) ? temp_arr[1] : temp_arr[2];
    }

    if (Verts[apex[1]].connect[0] == 3 || Verts[apex[2]].connect[0] == 3 || Verts[apex[3]].connect[0] == 3){
	//nop, we don't do this either
	return false;
    }



    //*******2) Destory one edge of the triangle and another one
    //there is 9 possibilities in near, we try them all :D

    int ip3 = ip[ip[0]];
    for (int i = 1; i <= ip[0]; i++){
	int j = (i == ip[0]) ? 1 : i + 1;

	int ip1 = ip[i];
	int ip2 = ip[j];
	int iapex = apex[i];
	int japex = apex[j];

	//we take the void_vertex as the avergae of two points
	//point one is the mid-point between ip1-ip2
	//point one is the mid-point between ip2-ip3
	void_vertex[0] = (Verts[ip1].x[0] + 2.0*Verts[ip2].x[0] + Verts[ip3].x[0]) / 4.0;
	void_vertex[1] = (Verts[ip1].x[1] + 2.0*Verts[ip2].x[1] + Verts[ip3].x[1]) / 4.0;
	void_vertex[2] = (Verts[ip1].x[2] + 2.0*Verts[ip2].x[2] + Verts[ip3].x[2]) / 4.0;


	//#1 first try to destroy the edges ip1-ip2 ip2-ip3
	if ((i == 1 || i == ip[0]) && Verts[ip2].connect[0] > 4){ //this prevents ip2 from being tri-valent
								  //(because we remove two verices connected to it and add one, effectively reduce valence by 1)

	    mynList[0] = 5;
	    mynList[1] = ip1;
	    mynList[2] = iapex;
	    mynList[3] = ip2;
	    mynList[4] = japex;
	    mynList[5] = ip3;

	    skipList[0] = 3;

	    if (inj){
		SetAttractorRepeller(Verts);
		AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
	    }

	    myConstraints.Reset(Verts, NULL);
	    myConstraints.Direction(NULL, mynList, Verts);

	    if (constraints.isDelaunay){
		//Delaunay constraints
		myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
	    }

	    if (myConstraints.OverlappingInSpheres()){
		if (constraints.isMinAngle || constraints.isMaxAngle){
		    //Min max angle
		    myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
		}

		if (constraints.isEdgeLen){
		    //Min edge length
		    myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
		}

		if (constraints.isSmooth){
		    //Smoothness
		    myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
		}


		StartActivePool(closedtSurfaceID, numSurfaceLayer);
		if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
		    //load edge to break
		    //ip1-ip2
		    //ip2-ip3
		    EdgeToBreak[0] = 2;
		    EdgeToBreak[1] = ip1;
		    EdgeToBreak[2] = ip2;
		    EdgeToBreak[3] = ip2;
		    EdgeToBreak[4] = ip3;
		    Inserter(numVert, Verts);
		    return true;
		}
	    }

	    if (inj){
		RevertAttractorRepeller(Verts);
	    }
	}



	//#2 second try to destroy the edges ip1-ip2 with another edge connected to on of the apex

	for (int id = 1; id <= 2; id++){
	    if (id == 2){
		std::swap(ip1, ip2);//trust me you need this.
	    }

	    if (Verts[ip1].connect[0] > 4){//this prevents ip1 from being tri-valent
					   //(because we remove to verices connected to it and add one, effectively reduce valence by 1)

					   //the two edges to destroy are ip1-ip2  and iapex1-ip1
		int extra_apex = FindCommonElement_SkipList(Verts[ip1].connect, Verts[iapex].connect, ip);

		if (i==2 && Verts[extra_apex].connect[0] > 3){//if it is 3-valent, then the void is ambigious
							      //the only solution would be creating 4-valent vertex which is equally bad for most cases
		    mynList[0] = 5;
		    mynList[1] = ip1;
		    mynList[2] = extra_apex;
		    mynList[3] = iapex;
		    mynList[4] = ip2;
		    mynList[5] = ip3;

		    //here we take the void_vertex as the avergae of two points
		    //point one is the mid-point between ip1-iapex
		    //point one is the mid-point between ip1-ip2
		    void_vertex[0] = (2.0*Verts[ip1].x[0] + Verts[iapex].x[0] + Verts[ip2].x[0]) / 4.0;
		    void_vertex[1] = (2.0*Verts[ip1].x[1] + Verts[iapex].x[1] + Verts[ip2].x[1]) / 4.0;
		    void_vertex[2] = (2.0*Verts[ip1].x[2] + Verts[iapex].x[2] + Verts[ip2].x[2]) / 4.0;

		    skipList[0] = 4;
		    skipList[skipList[0]] = iapex;

		    if (inj){
			SetAttractorRepeller(Verts);
			AttractorRepeller(ip, numVert, Verts, closedtSurfaceID, samplingBudget, numSurfaceLayer, false);
		    }

		    myConstraints.Reset(Verts, NULL);
		    myConstraints.Direction(NULL, mynList, Verts);

		    if (constraints.isDelaunay){
			//Delaunay constraints
			myConstraints.Delaunay(mynList, skipList, Verts, false, 1);
		    }

		    if (myConstraints.OverlappingInSpheres()){
			if (constraints.isMinAngle || constraints.isMaxAngle){
			    //Min max angle
			    myConstraints.MinMaxAngle(mynList, Verts, constraints.MinAngle, constraints.MaxAngle, void_vertex, constraints.isNonobtuse, 1);
			}

			if (constraints.isEdgeLen){
			    //Min edge length
			    myConstraints.MinEdgeLength(mynList, Verts, constraints.MinEdgeLength_sq);
			}

			if (constraints.isSmooth){
			    //Smoothness
			    myConstraints.Smoothness(mynList, skipList, Verts, void_vertex, constraints.dev);
			}


			StartActivePool(closedtSurfaceID, numSurfaceLayer);
			if (Sampler(Verts, mynList, samplingBudget, -1, -1, -1, &OpimizerFunc_CenterAngle)){
			    EdgeToBreak[0] = 2;
			    EdgeToBreak[1] = ip1;
			    EdgeToBreak[2] = ip2;
			    EdgeToBreak[3] = ip1;
			    EdgeToBreak[4] = iapex;
			    Inserter(numVert, Verts);
			    return true;
			}
		    }

		    if (inj){
			RevertAttractorRepeller(Verts);
		    }
		}
	    }
	}


	ip3 = ip2;//not ip1 because we swap it

    }



    return false;
}
void Operator::AttractorRepeller(int* ip, int&numVert, vert*Verts, int closedtSurfaceID, int samplingBudget, int numSurfaceLayer, bool isAttractor)
{
    //for each vertex in mynList, move it in eaither away (repeller) or closer (attractor)
    //in the direction of  void_vertex
    int prv_id = mynList[0];
    skipListNN[0] = 1;
    for (int i = 1; i <= mynList[0]; i++){
	int iq = mynList[i];
	//get the neighbour list of iq
	//this neighbour list is not a conected chain since iq is connected to ip's

	int nxt_id = (i == mynList[0]) ? 1 : i + 1;
	int start_id = GetIndex(mynList[prv_id], Verts[iq].connect);
	int end_id = GetIndex(mynList[nxt_id], Verts[iq].connect);

#ifdef DEBUGGING
	if (start_id < 0 || end_id < 0){
	    ErrWarnMessage(__LINE__, "Operator::Attractor:: can not grab the right index", 0);
	}
#endif

	mynListNN[0] = 0;
	if (isAttractor){
	    CycleList_SkipList(iq, Verts, ip, start_id, end_id, mynListNN);
	}
	else{
	    int before_prv = (prv_id == 1) ? mynList[mynList[0]] : mynList[prv_id - 1];
	    int after_nxt = (nxt_id == mynList[0]) ? mynList[1] : mynList[nxt_id + 1];

	    int move_away_vert;
	    if (GetIndex(before_prv, Verts[iq].connect) >= 0){
		move_away_vert = before_prv;
	    }
	    else if (GetIndex(after_nxt, Verts[iq].connect) >= 0){
		move_away_vert = after_nxt;
	    }
	    else{
		move_away_vert = mynList[nxt_id];
	    }

	    CycleList(iq, Verts, move_away_vert, start_id, end_id, mynListNN, 1);
	}
	if (mynListNN[0] < 0){
	    prv_id = i;
	    continue;
	}


	//build your geometric primitives to maintain the quality of triangles in mynListNN

	skipListNN[1] = iq;
	myConstraints.Reset(Verts, NULL);

	//myConstraints.Direction(skipListNN, mynListNN, Verts);
	myConstraints.Direction(skipListNN, Verts[iq].connect, Verts);

	if (constraints.isDelaunay){
	    if (!myConstraints.Delaunay(mynListNN, skipListNN, Verts, 1, 0)){
		prv_id = i;
		continue;
	    }

	    if (!myConstraints.DelaunayForcedNotConnected(mynList, iq, skipList, Verts)){
		prv_id = i;
		continue;
	    }

	    if (!myConstraints.OverlappingInSpheres()){
		prv_id = i;
		continue;
	    }
	}

	if (constraints.isMinAngle || constraints.isMaxAngle){
	    //Min max angle
	    myConstraints.MinMaxAngle(mynListNN, Verts, constraints.MinAngle, constraints.MaxAngle, Verts[iq].x, constraints.isNonobtuse, 0);
	}

	if (constraints.isEdgeLen){
	    //Min edge length
	    myConstraints.MinEdgeLength(mynListNN, Verts, constraints.MinEdgeLength_sq);
	}


	if (isAttractor){
	    myConstraints.AttractorInSphere(Verts, iq, void_vertex);
	}
	else{
	    //Hybird attractor/injection
	    double sum_angle = 0;
	    for (int kkk = 1; kkk < mynListNN[0]; kkk++){
		int n1 = mynListNN[kkk];
		int n2 = mynListNN[kkk + 1];
		sum_angle += AngleVectVect(Verts[n1].x[0] - Verts[iq].x[0], Verts[n1].x[1] - Verts[iq].x[1], Verts[n1].x[2] - Verts[iq].x[2],
			Verts[n2].x[0] - Verts[iq].x[0], Verts[n2].x[1] - Verts[iq].x[1], Verts[n2].x[2] - Verts[iq].x[2])*RadToDeg;
	    }
	    sum_angle = 360.0 - sum_angle;
	    if (sum_angle >= 0){
		if (sum_angle < 2.0*constraints.MinAngle){
		    myConstraints.AttractorInSphere(Verts, iq, void_vertex);
		}
		else if (sum_angle>2.0*constraints.MaxAngle){
		    myConstraints.RepellerExSphere(Verts, iq, void_vertex);
		}
		else{
		    //check to which it is closer (repeller or attracot)
		    double dif_min = sum_angle - 2.0*constraints.MinAngle;
		    double dif_max = 2.0*constraints.MaxAngle - sum_angle;
		    if (dif_max < dif_min){
			myConstraints.RepellerExSphere(Verts, iq, void_vertex);
		    }
		    else{
			myConstraints.AttractorInSphere(Verts, iq, void_vertex);
		    }
		}
	    }
	    else{
		//the default is to repel
		myConstraints.RepellerExSphere(Verts, iq, void_vertex);
	    }
	}

	//we construct a single plane such that the new point is never behind the the void vertex
	//it helps in case of large voids

	myConstraints.SinglePlane(Verts[iq].x[0] - void_vertex[0],
		Verts[iq].x[1] - void_vertex[1],
		Verts[iq].x[2] - void_vertex[2],
		void_vertex[0], void_vertex[1], void_vertex[2]);

	StartActivePool(closedtSurfaceID, numSurfaceLayer);
	if (Sampler(Verts, mynListNN, samplingBudget, mynList[nxt_id], iq, mynList[prv_id], (isAttractor) ? &OpimizerFunc_Closer : &OpimizerFunc_Further)){
	    Mover(iq, Verts);
	}

	prv_id = i;

    }
}
void Operator::SetAttractorRepeller(vert*Verts)
{
    //copy the vertices we gonna move for attracor or repeller
    for (int i = 1; i <= mynList[0]; i++){
	original_neighbout[i][0] = Verts[mynList[i]].x[0];
	original_neighbout[i][1] = Verts[mynList[i]].x[1];
	original_neighbout[i][2] = Verts[mynList[i]].x[2];
    }
}
void Operator::RevertAttractorRepeller(vert*Verts)
{
    //reset the vertices we may have moved for attractor or repeller
    //probably because we moved some of them and inejction or ejection failed
    for (int i = 1; i <= mynList[0]; i++){
	Verts[mynList[i]].x[0] = original_neighbout[i][0];
	Verts[mynList[i]].x[1] = original_neighbout[i][1];
	Verts[mynList[i]].x[2] = original_neighbout[i][2];
    }
}
void Operator::TriValent(int ip, vert*Verts, int*list, int skip1, int skip2)
{
    //store in list the trivalent nodes connected to ip
    //except if it skip1 or skip2
    for (int i = 1; i <= Verts[ip].connect[0]; i++){
	int iq = Verts[ip].connect[i];
	if (iq == skip1 || iq == skip2){ continue; }
	if (Verts[iq].connect[0] == 3){
	    list[++list[0]] = iq;
	}
    }
}
void Operator::GetListSkipVertex(int ip, vert*Verts, int skip1, int skip2, int*list)
{
    //store the connectivity list of ip in list with exception of skip1 and skip2
    list[0] = 0;
    for (int i = 1; i <= Verts[ip].connect[0]; i++){
	if (Verts[ip].connect[i] == skip1 || Verts[ip].connect[i] == skip2){ continue; }
	list[++list[0]] = Verts[ip].connect[i];
    }
}
void Operator::GetEdgeSortedNeighbourList(int ip1, int ip2, vert*Verts, int*list)
{
    //get the sorted list of neighbour connected to ip1 and ip2 such that
    //ip1 and ip2 are not tri-valent and they are not connected to shared tri-valent node

    //1) find the common node between ip1 and ip2 (should be two)
    int common[10];
    FindCommonElements(Verts[ip1].connect, Verts[ip2].connect, common);
    if (common[0] != 2){
	return;
    }

    int common1_id = GetIndex(common[1], Verts[ip1].connect);
    int common2_id = GetIndex(common[2], Verts[ip1].connect);

    //2) for ip1 connectivity list, start from common1 and keep advancing till you find common2
    //   for ip2 connectivity list, start from common2 and keep advancing till you find common1
    //   since sorting is not consistent, we could be advancing forward or backword in the connectivity list
    //   such that we are always moving away from the other ip

    list[0] = 0;

    CycleList(ip1, Verts, ip2, common1_id, common2_id, list, 0);
    if (list[0] < 0){ return; }

    common1_id = GetIndex(common[1], Verts[ip2].connect);
    common2_id = GetIndex(common[2], Verts[ip2].connect);

    CycleList(ip2, Verts, ip1, common2_id, common1_id, list, 0);


#ifdef DEBUGGING
    if (list[0] != Verts[ip1].connect[0] + Verts[ip2].connect[0] - 4){
	ErrWarnMessage(__LINE__, "Operator::GetEdgeSortedNeighbourList:: incorrect number of neighbours", 0);
    }
#endif

}
void Operator::GetFaceSortedNeighbourList(int ip1, int ip2, int ip3, vert*Verts, int*list)
{
    //get the sorted list of neighbour connected to ip1, ip2 and ip3 such that
    //ip1, ip2 and ip3 are not tri-valent and they are not connected to shared tri-valent node

    //1) find the common node between ip1 and ip2 (should be two)

    int ip1_ip2_sh = FindCommonElement_SkipList(Verts[ip1].connect, Verts[ip2].connect, ip3);
    int ip2_ip3_sh = FindCommonElement_SkipList(Verts[ip2].connect, Verts[ip3].connect, ip1);
    int ip3_ip1_sh = FindCommonElement_SkipList(Verts[ip3].connect, Verts[ip1].connect, ip2);

    if (ip1_ip2_sh < 0 || ip2_ip3_sh < 0 || ip3_ip1_sh < 0){
	list[0] = -1;
	return;
    }

    //2) for ip1 connectivity list, start from ip3_ip1_sh and keep advancing till you find ip1_ip2_sh
    //   for ip2 connectivity list, start from ip1_ip2_sh and keep advancing till you find ip2_ip3_sh
    //   for ip3 connectivity list, start from ip2_ip3_sh and keep advancing till you find ip3_ip1_sh
    //   since sorting is not consistent, we could be advancing forward or backword in the connectivity list
    //   such that when loop around ip1 we move away from ip3
    //             when loop around ip2 we move away from ip1
    //             when loop around ip3 we move away from ip2

    list[0] = 0;

    int sh12_id_1 = GetIndex(ip1_ip2_sh, Verts[ip1].connect);
    int sh31_id_1 = GetIndex(ip3_ip1_sh, Verts[ip1].connect);
    CycleList(ip1, Verts, ip3, sh31_id_1, sh12_id_1, list, 0);
    if (list[0] < 0){ return; }

    int sh12_id_2 = GetIndex(ip1_ip2_sh, Verts[ip2].connect);
    int sh23_id_2 = GetIndex(ip2_ip3_sh, Verts[ip2].connect);
    CycleList(ip2, Verts, ip1, sh12_id_2, sh23_id_2, list, 0);
    if (list[0] < 0){ return; }

    int sh23_id_3 = GetIndex(ip2_ip3_sh, Verts[ip3].connect);
    int sh31_id_3 = GetIndex(ip3_ip1_sh, Verts[ip3].connect);
    CycleList(ip3, Verts, ip2, sh23_id_3, sh31_id_3, list, 0);
    if (list[0] < 0){ return; }


#ifdef DEBUGGING
    //if (list[0] != Verts[ip1].connect[0] + Verts[ip2].connect[0] - 4){
    //	ErrWarnMessage(__LINE__, "Operator::GetFaceSortedNeighbourList:: incorrect number of neighbours", 0);
    //}
#endif

}
void Operator::CycleList_SkipList(int ip, vert*Verts, int*SkipList, int start_id, int end_id, int*list)
{
    //start_id and end_id are  id's in ip connectivity list
    //starting from start_id, move in ip connectivity list till we reach end_id
    //skip all vertices in SkipList which may or may not be in ip connectivity
    //store what you get in list
    //list should contain at the end a connected list of verices but only the start and end vertices are not connected

    int prv = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;

    prv = GetIndex(Verts[ip].connect[prv], SkipList);

    bool forward = (prv < 0);

    list[++list[0]] = Verts[ip].connect[start_id];
    int init_size = list[0];
    //predicate that we gonna move forward
    //if we meet any of ip's before reaching end_id, we than reset and go the other direction
    int my_start_id = start_id;
    while (true){
	my_start_id = (my_start_id == Verts[ip].connect[0]) ? 1 : my_start_id + 1;
	list[++list[0]] = Verts[ip].connect[my_start_id];
	if (GetIndex(Verts[ip].connect[my_start_id], SkipList) >= 0){ break; }
	if (my_start_id == end_id){ return; }
    }

    //if we reach here, then that means that moving forward was wrong
    list[0] = init_size;
    while (true){
	start_id = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;
	list[++list[0]] = Verts[ip].connect[start_id];
	if (start_id == end_id){ return; }
    }
}
void Operator::CycleList(int ip, vert*Verts, int iq, int start_id, int end_id, int*list, bool includeLastEle)
{
    //start_id and end_id are  id's in ip connectivity list
    //starting from start_id, move in ip connectivity list till we reach end_id
    //move such that we are moving away from iq

    int prv = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;
    int nxt = (start_id == Verts[ip].connect[0]) ? 1 : start_id + 1;

    if (Verts[ip].connect[prv] != iq && Verts[ip].connect[nxt] != iq){
	list[0] = -1;
	return;
	//ErrWarnMessage(__LINE__, "Operator::CycleList:: wrong id's", 0);
    }

    bool forward = (Verts[ip].connect[prv] == iq);

    list[++list[0]] = Verts[ip].connect[start_id];
    while (true){
	if (forward){
	    start_id = (start_id == Verts[ip].connect[0]) ? 1 : start_id + 1;
	}
	else{
	    start_id = (start_id == 1) ? Verts[ip].connect[0] : start_id - 1;
	}

	if (start_id == end_id){ break; }
	list[++list[0]] = Verts[ip].connect[start_id];
    }
    if (includeLastEle){
	list[++list[0]] = Verts[ip].connect[start_id];
    }
}
bool Operator::InspectFeasibleRegion(int*nList, vert*Verts)
{
    //here we inspect if there is feasible region
    //by inspect the void corners only (rather than relying on resampling to figure this out)

    //For the three operators, we form a void corner and then insert a vertex inside it
    //we check if a single vertex could possibly meet all the objectives or not
    //for example, if the void corners make an angle less than 2*min_angle
    //then any new vertex will split this angle to two angles such that one of them is less min_angle
    return true;

    if (constraints.isMinAngle || constraints.isMaxAngle){
	int n2 = nList[nList[0]];//previous
	for (int i = 1; i <= nList[0]; i++){
	    int ap = nList[i];//current
	    int n1 = (i == nList[0]) ? nList[1] : nList[i + 1]; //next
	    double xn, yn, zn;
	    //double angle = AngleVectVect(Verts[n1].x[0] - Verts[ap].x[0], Verts[n1].x[1] - Verts[ap].x[1], Verts[n1].x[2] - Verts[ap].x[2],
	    //	                         Verts[n2].x[0] - Verts[ap].x[0], Verts[n2].x[1] - Verts[ap].x[1], Verts[n2].x[2] - Verts[ap].x[2])*RadToDeg;//180.0 / PI;
	    double v1x = Verts[n1].x[0] - Verts[ap].x[0]; double v1y = Verts[n1].x[1] - Verts[ap].x[1]; double v1z = Verts[n1].x[2] - Verts[ap].x[2];
	    double v2x = Verts[n2].x[0] - Verts[ap].x[0]; double v2y = Verts[n2].x[1] - Verts[ap].x[1]; double v2z = Verts[n2].x[2] - Verts[ap].x[2];

	    Cross(v1x, v1y, v1z,
		    v2x, v2y, v2z,
		    xn, yn, zn);
	    NormalizeVector(xn, yn, zn);
	    double angle = Angle360Vectors(v1x, v1y, v1z, v2x, v2y, v2z, xn, yn, zn);
	    if (angle > 180.0){
		int sdfd = 4545;
	    }
	    if (constraints.isMinAngle && angle < 2.0*constraints.MinAngle){
		return false;
	    }
	    if (constraints.isMaxAngle && angle > 2.0*constraints.MaxAngle){
		return false;
	    }
	    n2 = ap;
	}
    }
    return true;

}

///************** Sampling
bool Operator::Sampler(vert*Verts, int*nList, int budget,
	int i_af, int i_rep, int i_bf,
	double(*Optimizer)(vert*Verts, int*nList, double*void_ver, double xx, double yy, double zz))
{

    //the resampling routine, called after constructing the geometric constraints
    //and initlize the active pool
    //it should spit out a single new vertex in case of success and return true
    //otherwise return false

    //here we should distingiush between optimizer sampling and regular sampling
    //optimizer seek to sample number of samples and pick the best to minimize the
    //returned value from Optimizer()
    //regular sampling just seek one sample to solve the problem
    bool att = i_rep >= 0;


    double currentBest = DBL_MAX;//we gonna minimize this
    if (budget > 1 && !att){
	//if optimizer, set the best to the current sample iff it meets all the constraints
	if (CheckNewVertex(Verts, nList, Verts[skipList[1]].x[0], Verts[skipList[1]].x[1], Verts[skipList[1]].x[2], att)){
	    currentBest = Optimizer(Verts, nList, void_vertex, Verts[skipList[1]].x[0], Verts[skipList[1]].x[1], Verts[skipList[1]].x[2]);

	    newVertex[0] = Verts[skipList[1]].x[0];
	    newVertex[1] = Verts[skipList[1]].x[1];
	    newVertex[2] = Verts[skipList[1]].x[2];
	}
    }

    int num_succ_candidate(0);

    for (int lf = 0; lf < 15; lf++){
	int attempt = int(0.8*num_active);

	//a) dart throwing

	for (int i = 0; i < attempt; i++){
	    int tri = int(myRandNum.RandNumGenerator()*double(num_active - 1));
	    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	    RetrieveCoordinates(lf, active_pool[tri], x1, y1, z1, x2, y2, z2, x3, y3, z3);
	    double x_new, y_new, z_new;
	    RandomPointInTri(x1, y1, z1, x2, y2, z2, x3, y3, z3, x_new, y_new, z_new);


	    if (CheckNewVertex(Verts, nList, x_new, y_new, z_new, att)){
		num_succ_candidate++;
		if (budget == 1){
		    newVertex[0] = x_new;
		    newVertex[1] = y_new;
		    newVertex[2] = z_new;
		    return true;
		}
		else{
		    //apply optimizer
		    double myVal = Optimizer(Verts, nList, void_vertex, x_new, y_new, z_new);
		    if (abs(currentBest - DBL_MAX) < _tol){
			//if this is the first candidate (then take it)
			currentBest = myVal;
			newVertex[0] = x_new;
			newVertex[1] = y_new;
			newVertex[2] = z_new;
		    }
		    else{

			if (myVal < currentBest){
			    currentBest = myVal;
			    newVertex[0] = x_new;
			    newVertex[1] = y_new;
			    newVertex[2] = z_new;
			}
		    }
		    if (num_succ_candidate >= budget){
			return true;
		    }
		}
	    }
	}

	//b) refinement
	int tmp_num_active = 0;
	for (int iactive = 0; iactive < num_active; iactive++){
	    double x1, y1, z1, x2, y2, z2, x3, y3, z3;
	    RetrieveCoordinates(lf, active_pool[iactive], x1, y1, z1, x2, y2, z2, x3, y3, z3);

	    int C = 0;
	    _tar[C][0] = x1;          _tar[C][1] = y1;          _tar[C][2] = z1;
	    _tar[C][3] = (x1 + x2) / 2.0; _tar[C][4] = (y1 + y2) / 2.0; _tar[C][5] = (z1 + z2) / 2.0;
	    _tar[C][6] = (x1 + x3) / 2.0; _tar[C][7] = (y1 + y3) / 2.0; _tar[C][8] = (z1 + z3) / 2.0;
	    C = 1;
	    _tar[C][0] = (x1 + x2) / 2.0; _tar[C][1] = (y1 + y2) / 2.0; _tar[C][2] = (z1 + z2) / 2.0;
	    _tar[C][3] = x2;          _tar[C][4] = y2;          _tar[C][5] = z2;
	    _tar[C][6] = (x2 + x3) / 2.0; _tar[C][7] = (y2 + y3) / 2.0; _tar[C][8] = (z2 + z3) / 2.0;
	    C = 2;
	    _tar[C][0] = (x1 + x3) / 2.0; _tar[C][1] = (y1 + y3) / 2.0; _tar[C][2] = (z1 + z3) / 2.0;
	    _tar[C][3] = (x2 + x3) / 2.0; _tar[C][4] = (y2 + y3) / 2.0; _tar[C][5] = (z2 + z3) / 2.0;
	    _tar[C][6] = x3;          _tar[C][7] = y3;          _tar[C][8] = z3;
	    C = 3;
	    _tar[C][0] = (x1 + x3) / 2.0; _tar[C][1] = (y1 + y3) / 2.0; _tar[C][2] = (z1 + z3) / 2.0;
	    _tar[C][3] = (x1 + x2) / 2.0; _tar[C][4] = (y1 + y2) / 2.0; _tar[C][5] = (z1 + z2) / 2.0;
	    _tar[C][6] = (x2 + x3) / 2.0; _tar[C][7] = (y2 + y3) / 2.0; _tar[C][8] = (z2 + z3) / 2.0;

	    for (C = 0; C < 4; C++){


		if (CheckRefinement(Verts, nList, _tar[C])){

		    tmp_active_pool[tmp_num_active][0] = active_pool[iactive][0];
		    tmp_active_pool[tmp_num_active][1] = active_pool[iactive][1];
		    tmp_active_pool[tmp_num_active][2] = active_pool[iactive][2];

		    if (lf>15){
			size_t shifted = C << 2 * (lf - 16);
			tmp_active_pool[tmp_num_active][2] = tmp_active_pool[tmp_num_active][2] | shifted;
		    }
		    else{
			size_t shifted = C << (2 * (lf));
			tmp_active_pool[tmp_num_active][1] = tmp_active_pool[tmp_num_active][1] | shifted;
		    }
		    tmp_num_active += 1;
		    if (tmp_num_active + 2 > MAXPOOL){
			//cout<<"Too much refining in Resampling()"<<endl;
			if (num_succ_candidate > 0){ return true; }
			return false;

		    }
		}
	    }
	}

	//swap
	if (tmp_num_active == 0){
	    if (num_succ_candidate > 0){ return true; }
	    return false;
	}
	else{
	    num_active = tmp_num_active;
	    for (int iactive = 0; iactive < num_active; iactive++){
		active_pool[iactive][0] = tmp_active_pool[iactive][0];
		active_pool[iactive][1] = tmp_active_pool[iactive][1];
		active_pool[iactive][2] = tmp_active_pool[iactive][2];
	    }
	}
    }

    return false;//why did we get here !!!

}
void Operator::StartActivePool(int id, int maxDepth)
{
    //TODO pass me the kd-tree of the original surface and let me pick the closest the surface myself

    //Start the active pool by the fan triangle of id
    //then propagate to include maxDepth of surrounding triangles of above fan triangle

    //for relocation, it is better to let maxDepth to be 1
    //for ejection, it is better to get larger patch tho 2 would be enough

    num_active = 0;
    next_layer[0] = 1;
    next_layer[1] = id;
    int start = 1;
    for (int lay = 0; lay < maxDepth; lay++){
	int end = next_layer[0];
	for (int i = start; i <=end; i++){
	    int myID = next_layer[i];
	    for (int j = 1; j <= Vert_org[myID].connect[0]; j++){
		int candid = Vert_org[myID].connect[j];
		if (GetIndex(candid, next_layer) < 0){
		    next_layer[++next_layer[0]] = candid;
#ifdef DEBUGGING
		    if (next_layer[0] + 1 > MAXPOOL){
			ErrWarnMessage(__LINE__, "Operator::StartActivePool:: increase length of next_layer", 0);
		    }
#endif
		}
	    }
	}
	start = end + 1;
    }



    for (int i = 1; i <= next_layer[0]; i++){
	GetFanTriangles(next_layer[i]);
    }

}
void Operator::StartActivePool(int closedtSurfaceID, double*myVert)
{
    //initialize actove pool but triangle fans of some surface vertices
    //find the surface trianle such that the projection of myVert is min
    //get the fan triangle of the three vertices of this triangle

    //myVert is the coordinates of the center to be projected
    //for relocation, it is the vertex to be reolcated
    //for ejection/injection, it the mid point of the two points to be ejected, or center of triangle to be ejected
    num_active = 0;

    int closedtSurfaceID1, closedtSurfaceID2;
    int n2 = Vert_org[closedtSurfaceID].connect[Vert_org[closedtSurfaceID].connect[0]];
    double closestDist = DBL_MAX;
    for (int i = 1; i <= Vert_org[closedtSurfaceID].connect[0]; i++){
	int n1 = Vert_org[closedtSurfaceID].connect[i];
	double x_projected, y_projected, z_projected;
	double dist = PointTriangleDistance(Vert_org[closedtSurfaceID].x[0], Vert_org[closedtSurfaceID].x[1], Vert_org[closedtSurfaceID].x[2],
		Vert_org[n1].x[0], Vert_org[n1].x[1], Vert_org[n1].x[2],
		Vert_org[n2].x[0], Vert_org[n2].x[1], Vert_org[n2].x[2],
		myVert[0], myVert[1], myVert[2],
		x_projected, y_projected, z_projected);
	if (dist < closestDist){
	    closestDist = dist;
	    closedtSurfaceID1 = n1;
	    closedtSurfaceID2 = n2;
	}
	n2 = n1;
    }


    //get the fan triangles of closedtSurfaceID, closedtSurfaceID1 and closedtSurfaceID2
    //carful of duplication
    GetFanTriangles(closedtSurfaceID);
    GetFanTriangles(closedtSurfaceID1);
    GetFanTriangles(closedtSurfaceID2);

}
void Operator::GetFanTriangles(int id)
{
    //get the fan triangles in id and store them in active_pool
    //make sure there is no duplication
    int n2 = Vert_org[id].connect[Vert_org[id].connect[0]];

    for (int i = 1; i <= Vert_org[id].connect[0]; i++){

	int n1 = Vert_org[id].connect[i];

	//id,n1,n2
	//check duplication
	bool dup = false;
	for (int j = 0; j < num_active; j++){
	    if (IsDuplicated(n2, n1, id, j)){
		dup = true;
		break;
	    }
	}

	if (!dup){
	    tri_pool[num_active][0] = id;
	    tri_pool[num_active][1] = n1;
	    tri_pool[num_active][2] = n2;

	    active_pool[num_active][0] = num_active;
	    active_pool[num_active][1] = 0;
	    active_pool[num_active][2] = 0;
	    num_active++;

	    if (num_active + 2 >= MAXPOOL){ return; }

	}

	n2 = n1;
    }
}
bool Operator::IsDuplicated(int ip1, int ip2, int ip3, int t2)
{
    return (ip1 == tri_pool[t2][0] || ip1 == tri_pool[t2][1] || ip1 == tri_pool[t2][2]) &&
	(ip2 == tri_pool[t2][0] || ip2 == tri_pool[t2][1] || ip2 == tri_pool[t2][2]) &&
	(ip3 == tri_pool[t2][0] || ip3 == tri_pool[t2][1] || ip3 == tri_pool[t2][2]);

}
void Operator::RandomPointInTri(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double&x, double&y, double&z)
{
    double r1, r2;

    r1 = myRandNum.RandNumGenerator();
    r2 = myRandNum.RandNumGenerator();

    x = (1 - sqrt(r1))*x1 + (sqrt(r1)*(1 - r2))*x2 + (r2*sqrt(r1))*x3;
    y = (1 - sqrt(r1))*y1 + (sqrt(r1)*(1 - r2))*y2 + (r2*sqrt(r1))*y3;
    z = (1 - sqrt(r1))*z1 + (sqrt(r1)*(1 - r2))*z2 + (r2*sqrt(r1))*z3;

}
bool Operator::CheckRefinement(vert*Verts, int*nList, double*tri)
{


    if (!myConstraints.InsideFeasibleRegion_Triangle(tri)){ return false; }

    /*if (constraints.isMinAngle || constraints.isMaxAngle){ //this is wrong
							     //TODO include the attractor here
							     int n2 = nList[nList[0]];
							     for (int i = 1; i <= nList[0]; i++){
							     int n1 = nList[i];
							     double 	angle1 = (AngleVectVect(Verts[n1].x[0] - tri[0], Verts[n1].x[1] - tri[1], Verts[n1].x[2] - tri[2],
							     Verts[n2].x[0] - tri[0], Verts[n2].x[1] - tri[1], Verts[n2].x[2] - tri[2]))*RadToDeg;

							     double angle2 = (AngleVectVect(Verts[n1].x[0] - tri[3], Verts[n1].x[1] - tri[4], Verts[n1].x[2] - tri[5],
							     Verts[n2].x[0] - tri[3], Verts[n2].x[1] - tri[4], Verts[n2].x[2] - tri[5]))*RadToDeg;

							     double angle3 = (AngleVectVect(Verts[n1].x[0] - tri[6], Verts[n1].x[1] - tri[7], Verts[n1].x[2] - tri[8],
							     Verts[n2].x[0] - tri[6], Verts[n2].x[1] - tri[7], Verts[n2].x[2] - tri[8]))*RadToDeg;
							     if ((angle1 < constraints.MinAngle + _tol &&
							     angle2 < constraints.MinAngle + _tol &&
							     angle3 < constraints.MinAngle + _tol)||
							     (angle1 > constraints.MaxAngle + _tol &&
							     angle2 > constraints.MaxAngle + _tol &&
							     angle3 > constraints.MaxAngle + _tol)){
							     return false;
							     }
							     n2 = n1;
							     }
							     }*/
    return true;

}

bool Operator::CheckNewVertex(vert*Verts, int*nList, double x_new, double  y_new, double  z_new, bool att)
{
    //chech if a new vertex is an acceptable vertex (outside all exclusion regions,
    //and inside all inclusion regions)
    if (!myConstraints.InsideFeasibleRegion_Vertex(x_new, y_new, z_new)){ return false; }

    if (constraints.isMinAngle || constraints.isMaxAngle){
	//check the apex angle

	if (att){
	    double accum = 0;
	    for (int i = 1; i < nList[0]; i++){
		int n2 = nList[i];
		int n1 = nList[i + 1];
		double myAngle = AngleVectVect(Verts[n1].x[0] - x_new, Verts[n1].x[1] - y_new, Verts[n1].x[2] - z_new,
			Verts[n2].x[0] - x_new, Verts[n2].x[1] - y_new, Verts[n2].x[2] - z_new)*RadToDeg;

		if (myAngle < constraints.MinAngle + _tol){ return false; }
		if (myAngle + _tol > constraints.MaxAngle){ return false; }
		accum += myAngle;
	    }
	}
	else{
	    int n2 = nList[nList[0]];
	    for (int i = 1; i <= nList[0]; i++){
		int n1 = nList[i];

		double myAngle = AngleVectVect(Verts[n1].x[0] - x_new, Verts[n1].x[1] - y_new, Verts[n1].x[2] - z_new,
			Verts[n2].x[0] - x_new, Verts[n2].x[1] - y_new, Verts[n2].x[2] - z_new)*RadToDeg;

		if (myAngle < constraints.MinAngle + _tol || myAngle + _tol > constraints.MaxAngle){ return false; }
		n2 = n1;
	    }
	}
    }

    if (constraints.isSmooth){
	//TODO add smoothness to attractor

	//check smoothnes
	int n2 = nList[nList[0]];
	for (int i = 1; i <= nList[0]; i++){
	    int n1 = nList[i];
	    int n3 = (i == nList[0]) ? nList[1] : nList[i + 1];
	    double angle = TriTriNormalAngle3(x_new, y_new, z_new, Verts[n1].x, Verts[n2].x, Verts[n3].x);
	    //angle = 180 - angle;
	    //if (angle < 180 - constraints.dev) {
	    if (angle > constraints.dev){
		return false;
	    }
	    n2 = n1;
	}
    }
    return true;
}

void Operator::RetrieveCoordinates(int lf, int* cell, double&x0, double&y0, double&z0, double&x1, double&y1, double&z1, double&x2, double&y2, double&z2)
{
    int po = tri_pool[cell[0]][0];
    int p1 = tri_pool[cell[0]][1];
    int p2 = tri_pool[cell[0]][2];

    x0 = Vert_org[po].x[0];
    y0 = Vert_org[po].x[1];
    z0 = Vert_org[po].x[2];
    x1 = Vert_org[p1].x[0];
    y1 = Vert_org[p1].x[1];
    z1 = Vert_org[p1].x[2];
    x2 = Vert_org[p2].x[0];
    y2 = Vert_org[p2].x[1];
    z2 = Vert_org[p2].x[2];

    size_t binary = cell[1], c_num;//c_num is the ID of refiend triangle
				   //binary is the avariable holding the value of the cell[1 or 2] which will be "bit" masked to get the refined ID
    int shift = 0;
    double a0, b0, c0, a1, b1, c1, a2, b2, c2;
    for (int i = 0; i < lf; i++){
	//get refined tri ID
	c_num = binary& ((int(pow(2.0, (shift + 1.0))) + 1) | (int(pow(2.0, shift)) + 1));
	c_num = c_num >> shift;
	//cout<<c_num<<"\n";
	//switch and get new coordinates of C
	switch (c_num)
	{
	    case 0:

		x1 = (x0 + x1) / 2.0;
		y1 = (y0 + y1) / 2.0;
		z1 = (z0 + z1) / 2.0;

		x2 = (x0 + x2) / 2.0;
		y2 = (y0 + y2) / 2.0;
		z2 = (z0 + z2) / 2.0;
		break;
	    case 1:

		x0 = (x0 + x1) / 2.0; y0 = (y0 + y1) / 2.0; z0 = (z0 + z1) / 2.0;

		x2 = (x1 + x2) / 2.0; y2 = (y1 + y2) / 2.0; z2 = (z1 + z2) / 2.0;
		break;
	    case 2:

		x0 = (x0 + x2) / 2.0; y0 = (y0 + y2) / 2.0; z0 = (z0 + z2) / 2.0;

		x1 = (x1 + x2) / 2.0; y1 = (y1 + y2) / 2.0; z1 = (z1 + z2) / 2.0;
		break;
	    case 3://in c3 I pick Po as the (right side direction point) and then counter clockwise

		a0 = (x0 + x2) / 2.0; b0 = (y0 + y2) / 2.0; c0 = (z0 + z2) / 2.0;

		a1 = (x0 + x1) / 2.0; b1 = (y0 + y1) / 2.0; c1 = (z0 + z1) / 2.0;

		a2 = (x1 + x2) / 2.0; b2 = (y1 + y2) / 2.0; c2 = (z1 + z2) / 2.0;

		x0 = a0; y0 = b0; z0 = c0; x1 = a1; y1 = b1; z1 = c1; x2 = a2; y2 = b2; z2 = c2;
		break;
	    default:
		fprintf(stderr, "\nError.\n");
	}

	if (i == 15){//switching to the second cell after 16 itaration of refinement
	    binary = cell[2];
	    shift = 0;
	}
	shift = shift + 2;
    }

}


///************** Optimizer Functions (for sampling)
double OpimizerFunc_CenterAngle(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz)
{
    //return the different between myAngle and perfect_angle
    //perfect_angle = 360.0 / nList[0];
    //myAngle is the min apex angle

    double perfect_angle = 360.0 / nList[0];

    double myAngle = 360.0;
    int n2 = nList[nList[0]];
    for (int i = 1; i <= nList[0]; i++){
	int n1 = nList[i];
	double ang = AngleVectVect(Verts[n1].x[0] - xx, Verts[n1].x[1] - yy, Verts[n1].x[2] - zz,
		Verts[n2].x[0] - xx, Verts[n2].x[1] - yy, Verts[n2].x[2] - zz)*RadToDeg;
	myAngle = std::min(myAngle, ang);
	n2 = n1;
    }
    return abs(myAngle - perfect_angle);
}
double OpimizerFunc_SideAngle(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz)
{
    //get the angle between (xx,yy,zz), nList[1] and nList[nList[0]]
    //if it supposed to be less than 180.0
    //otherwise this method might returen wrong angle

    //return the difference between (2minAngle - ang) + (ang - 2maxAngle)
    //we take the difference this way because we try to minimize the returned value from this function

    int n1 = nList[1];
    int n2 = nList[nList[0]];
    double ang = AngleVectVect(Verts[n1].x[0] - xx, Verts[n1].x[1] - yy, Verts[n1].x[2] - zz,
	    Verts[n2].x[0] - xx, Verts[n2].x[1] - yy, Verts[n2].x[2] - zz)*RadToDeg;


    double diff1 = 2.0*minAngle - ang;
    double diff2 = ang - 2.0*maxAngle;


    return diff2 + diff2;
}
double OpimizerFunc_SideAngleHybird(vert*Verts, int*nList, double minAngle, double maxAngle, double xx, double yy, double zz)
{
    //TODO (smarter way for optimizing side angles)
    //use this for hybird attractor/repeller

    //nList is discounted list of neighbours around (xx,yy,zz)
    //get sum of all apex angle at (xx,yy,zz)
    //subtract this from 360 --> sum_angle
    //depending min(diff_min, diff_max)
    //where diff_min = sum_angle - 2.0*minAngle
    //diff_max = 2.0*maxAngle - sum_angle
    double sum_angle = 0;
    for (int i = 1; i < nList[0]; i++){
	int n1 = nList[i];
	int n2 = nList[i + 1];
	sum_angle += AngleVectVect(Verts[n1].x[0] - xx, Verts[n1].x[1] - yy, Verts[n1].x[2] - zz,
		Verts[n2].x[0] - xx, Verts[n2].x[1] - yy, Verts[n2].x[2] - zz)*RadToDeg;
    }
    sum_angle = 360.0 - sum_angle;

    double diff_min = -(sum_angle - 2.0*minAngle);
    double diff_max = -(2.0*maxAngle - sum_angle);

    return std::min(diff_min, diff_max);

}
double OpimizerFunc_Closer(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz)
{
    //return the distance to void_vertex
    //since we optimizer minimize the distance, we return it as it is

    double dist = Dist(void_verx[0], void_verx[1], void_verx[2], xx, yy, zz);
    return dist;
}
double OpimizerFunc_Further(vert*Verts, int*nList, double*void_verx, double xx, double yy, double zz)
{
    //return the distance to void_vertex
    //since we optimizer minimize the distance, we return negative the distance

    double dist = Dist(void_verx[0], void_verx[1], void_verx[2], xx, yy, zz);
    return -dist;
}



///************** Update Mesh
void Operator::Mover(int ip, vert*Verts)
{
    //move ip to newVertex (do any other neccessary toplogical changes)
    Verts[ip].x[0] = newVertex[0];
    Verts[ip].x[1] = newVertex[1];
    Verts[ip].x[2] = newVertex[2];
}
void Operator::EdgeCollapse(int ip_low, int ip_high, int&numVert, vert*Verts)
{
    //collapse the edge between ip_low and ip_high to a point (newVertex)
    //remove ip_high and then move ip_low (the order of the operation is i m p o r t a nt)
    if (ip_low > ip_high){ std::swap(ip_high, ip_low); }

    //1) remove ip_high
    RemoveVertex(ip_high, numVert, Verts);

    //2) move ip_low
    Verts[ip_low].x[0] = newVertex[0];
    Verts[ip_low].x[1] = newVertex[1];
    Verts[ip_low].x[2] = newVertex[2];
    //remove ip_low from any connectivity list of other vertices
    for (int i = 1; i <= Verts[ip_low].connect[0]; i++){
	int iq = Verts[ip_low].connect[i];
	RemoveNodeFromList(Verts[iq].connect, ip_low);
    }

    //update the connectivity of ip_low to be mynList
    for (int i = 0; i <= mynList[0]; i++){
	Verts[ip_low].connect[i] = mynList[i];
    }

    //appropriately, add ip_low to the connectivity list of vertices on mynList
    int prv = mynList[mynList[0]];//prv
    for (int i = 1; i <= mynList[0]; i++){
	int nxt = (i == mynList[0]) ? mynList[1] : mynList[i + 1]; //next
	AddEntrySortedList(ip_low, Verts[mynList[i]].connect, prv, nxt);
	prv = mynList[i];
    }


}
void Operator::FaceCollapse(int ip_low, int ip_mid, int ip_high, int&numVert, vert*Verts)
{
    //collapse the triangle ip_low-ip_mid-ip_high
    //remove ip_high and then remove ip_mid and then move ip_low (the order of the operation is i m p o r t a nt)

    int my_ip_high, my_ip_mid, my_ip_low;
    my_ip_high = std::max(ip_low, ip_mid);
    my_ip_high = std::max(my_ip_high, ip_high);

    my_ip_low = std::min(ip_low, ip_mid);
    my_ip_low = std::min(my_ip_low, ip_high);

    my_ip_mid = (ip_mid > my_ip_low && ip_mid<my_ip_high) ? ip_mid : ((ip_low>my_ip_low && ip_low < my_ip_high) ? ip_low : ip_high);

    ip_low = my_ip_low; ip_high = my_ip_high; ip_mid = my_ip_mid;

    //1) remove ip_high, ip_mid
    RemoveVertex(ip_high, numVert, Verts);
    RemoveVertex(ip_mid, numVert, Verts);

    //2) move ip_low
    Verts[ip_low].x[0] = newVertex[0];
    Verts[ip_low].x[1] = newVertex[1];
    Verts[ip_low].x[2] = newVertex[2];
    //remove ip_low from any connectivity list of other vertices
    for (int i = 1; i <= Verts[ip_low].connect[0]; i++){
	int iq = Verts[ip_low].connect[i];
	RemoveNodeFromList(Verts[iq].connect, ip_low);
    }

    //update the connectivity of ip_low to be mynList
    for (int i = 0; i <= mynList[0]; i++){
	Verts[ip_low].connect[i] = mynList[i];
    }

    //appropriately, add ip_low to the connectivity list of vertices on mynList
    int prv = mynList[mynList[0]];//prv
    for (int i = 1; i <= mynList[0]; i++){
	int nxt = (i == mynList[0]) ? mynList[1] : mynList[i + 1]; //next
	AddEntrySortedList(ip_low, Verts[mynList[i]].connect, prv, nxt);
	prv = mynList[i];
    }

}
void Operator::RemoveVertex(int ip, int&numVert, vert*Verts)
{
    //remove ip from mesh data strcuture

    //1) remove ip from any connectivity list of other vertices
    for (int i = 1; i <= Verts[ip].connect[0]; i++){
	int iq = Verts[ip].connect[i];
	RemoveNodeFromList(Verts[iq].connect, ip);
    }

    //2)grab the last mesh vertex and promote it to be ip
    //if ip was already the last mesh vertex, then never mind this :)
    numVert--;
    if (ip!=numVert){
	for (int i = 0; i <= Verts[numVert].connect[0]; i++){
	    Verts[ip].connect[i] = Verts[numVert].connect[i];
	}
	Verts[ip].x[0] = Verts[numVert].x[0];
	Verts[ip].x[1] = Verts[numVert].x[1];
	Verts[ip].x[2] = Verts[numVert].x[2];

	//3) any vertex connected to numVert should know it new name -> ip
	for (int i = 1; i <= Verts[numVert].connect[0]; i++){
	    int iq = Verts[numVert].connect[i];
	    int id = GetIndex(numVert, Verts[iq].connect);
#ifdef DEBUGGING
	    if (id < 0){
		ErrWarnMessage(__LINE__, "Operator::RemoveVertex error(1)", 0);
	    }
#endif
	    Verts[iq].connect[id] = ip;
	}
    }

    //numVert could be also there in mynList (which will be used in update another mesh connectivity in case of edge collapse)
    //so, we gotta update this one too
    int id = GetIndex(numVert, mynList);
    if (id >= 0){
	mynList[id] = ip;
    }


}
void Operator::Inserter(int&numVert, vert*Verts)
{
    //break edges in EdgeToBreak,
    //insert a newVertex into the data structure and connect it to mynList,
    //update mynList too,

    for (int i = 1; i <= EdgeToBreak[0]; i++){
	int ip = EdgeToBreak[i * 2 - 1];
	int iq = EdgeToBreak[i * 2];
	RemoveNodeFromList(Verts[iq].connect, ip);
	RemoveNodeFromList(Verts[ip].connect, iq);
    }

    Verts[numVert].x[0] = newVertex[0];
    Verts[numVert].x[1] = newVertex[1];
    Verts[numVert].x[2] = newVertex[2];

    for (int i = 0; i <= mynList[0]; i++){
	Verts[numVert].connect[i] = mynList[i];
    }


    //appropriately, add numVert  to the connectivity list of vertices on mynList
    int prv = mynList[mynList[0]];//prv
    for (int i = 1; i <= mynList[0]; i++){
	int nxt = (i == mynList[0]) ? mynList[1] : mynList[i + 1]; //next
	AddEntrySortedList(numVert, Verts[mynList[i]].connect, prv, nxt);
	prv = mynList[i];
    }

    numVert++;

}


///************** Debug
void Operator::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{
    //mess_id =0 for error (exit)
    //otherwise, it is a warning (pause)

    if (mess_id == 0){
	fprintf(stderr, "\nError::line(%d)-->>%s", lineNum, message.c_str());
	system("pause");
    }
    else{
	fprintf(stderr, "\nWarning::line(%d)-->>%s\n", lineNum, message.c_str());
	system("pause");
    }
}
void Operator::DrawActivePool(int lf)
{
    std::cout << "\n I AM DRAWING in Operator::DrawActivePool()" << std::endl;

    std::stringstream fname;
    fname << "debug_out/active_pool.obj";
    std::fstream file(fname.str().c_str(), std::ios::out);
    double x1, y1, z1, x2, y2, z2, x3, y3, z3;


    for (int V = 0; V<num_active; V++){
	RetrieveCoordinates(lf, active_pool[V], x1, y1, z1, x2, y2, z2, x3, y3, z3);
	file << "v " << x1 << " " << y1 << " " << z1 << std::endl;
	file << "v " << x2 << " " << y2 << " " << z2 << std::endl;
	file << "v " << x3 << " " << y3 << " " << z3 << std::endl;
    }
    for (int V = 1; V <= 3 * num_active; V += 3){
	file << "f " << V << " " << V + 1 << " " << V + 2 << std::endl;
    }
}
void Operator::DrawnList(vert*Verts)
{
    for (int i = 1; i <= mynList[0]; i++){
	//DrawOneSphere(mynList[i], Verts[mynList[i]].x[0], Verts[mynList[i]].x[1], Verts[mynList[i]].x[2], 0.001, 3);
    }

}

void CMeshImp::ErrWarnMessage(size_t lineNum, std::string message, size_t mess_id)
{
    //mess_id =0 for error (exit)
    //otherwise, it is a warning (pause)

    if (mess_id == 0){
	fprintf(stderr, "\nError::line(%d)-->>%s", lineNum, message.c_str());
	exit(1);
    }
    else{
	fprintf(stderr, "\nWarning::line(%d)-->>%s\n", lineNum, message.c_str());
	system("pause");
    }
}

//** Containers **//
CMeshImp::CMeshImp(int numVert, double **Verts, int numTri, int**Tris) : numVert_org(numVert), numTri_org(numTri)
{

    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "*************************** Initialize Data Structure ***************************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");


    //initiate containers
    MaxNumVert = size_t(3.5*numVert); //TODO only set the 1.5 factor when the injection operators are used
    Vert_org = new vert[MaxNumVert];
    Vert_imp = new vert[MaxNumVert];
    numVert_imp = numVert_org;

    myBoundingBox.xmax = myBoundingBox.ymax = myBoundingBox.zmax = DBL_MIN;
    myBoundingBox.xmin = myBoundingBox.ymin = myBoundingBox.zmin = DBL_MAX;


    //1) store vertices coordinates and set bounding  box
    for (int i = 0; i < numVert_org; i++){
	for (int j = 0; j < 3; j++){
	    Vert_imp[i].x[j] = Vert_org[i].x[j] = Verts[i][j];
	}
	Vert_imp[i].connect[0] = Vert_org[i].connect[0] = 0;

	myBoundingBox.xmax = std::max(myBoundingBox.xmax, Verts[i][0]);
	myBoundingBox.xmin = std::min(myBoundingBox.xmin, Verts[i][0]);
	myBoundingBox.ymax = std::max(myBoundingBox.ymax, Verts[i][1]);
	myBoundingBox.ymin = std::min(myBoundingBox.ymin, Verts[i][1]);
	myBoundingBox.zmax = std::max(myBoundingBox.zmax, Verts[i][2]);
	myBoundingBox.zmin = std::min(myBoundingBox.zmin, Verts[i][2]);
    }

    myBoundingBox.lx = myBoundingBox.xmax - myBoundingBox.xmin;
    myBoundingBox.ly = myBoundingBox.ymax - myBoundingBox.ymin;
    myBoundingBox.lz = myBoundingBox.zmax - myBoundingBox.zmin;


    //2) store  tessellation  and find neighbour triangles
    Tri_org = new tri[numTri_org];
    for (int i = 0; i < numTri_org; i++){
	for (int j = 0; j < 3; j++){
	    Tri_org[i].id[j] = Tris[i][j];
	}
    }
    int**vertex_triangles = NULL; //for each vertex i, what are the triangles i is a head on  (this array is populated in FindNeighbourTriangles)
    FindNeighbourTriangles(vertex_triangles);

    //3) scale inside a unit box
    ScaleInputMesh();

    //4) Build KdTree for the input surface
    surfaceKdTree.BuildTree(numVert, Vert_org, 3);

    //5) get connectivity for Vert_imp/Vert_org (fan triangle)
    GetFanTriangle();

    //6) check for isolated vertices
    for (int i = 0; i < numVert_org; i++){
	if (Vert_org[i].connect[0] == 0){
	    ErrWarnMessage(__LINE__, "CMeshImp::CMeshImp:: isolated vertex id=" + std::to_string(i), 1);
	}
    }

    //7) sort the connectivity
    InitialSortAroundVertex(vertex_triangles);
    delete[]vertex_triangles;


    //8)Get Statistics
    DisplayStat(numVert_org, Vert_org, 1);
    //WriteStatToFile("input_stat.txt", numVert_org, Vert_org, 1);

    //9)Draw input mesh
    //GetMeshOBJ("input.obj", 1, "input_obtuse.obj");
}
CMeshImp::~CMeshImp(){};
void CMeshImp::FindNeighbourTriangles(int**&myTriangles)
{
    myTriangles = new int*[numVert_org]; //for each vertex i, temporary store the triangles where i is a head on
    int common[20];
    for (int i = 0; i < numVert_org; i++){
	myTriangles[i] = new int[20];
	myTriangles[i][0] = 0;
    }

    //loop over triangles and populate myTriangles
    for (int t = 0; t < numTri_org; t++){
	int id0(Tri_org[t].id[0]), id1(Tri_org[t].id[1]), id2(Tri_org[t].id[2]);
	myTriangles[id0][++myTriangles[id0][0]] = t;
	myTriangles[id1][++myTriangles[id1][0]] = t;
	myTriangles[id2][++myTriangles[id2][0]] = t;

	if (myTriangles[id0][0] >= 20 || myTriangles[id1][0] >= 20 || myTriangles[id2][0] >= 20){
	    ErrWarnMessage(__LINE__, "CMeshImp::FindNeighbourTriangles:: a vertex is shared within more than 20 triangles!!", 0);
	}
    }

    for (int t = 0; t < numTri_org; t++){
	//for each triangle t
	//find the other triangle shared between two of its heads
	//update the t's neighbour list

	for (int i = 0; i < 3; i++){
	    int j = (i == 2) ? 0 : i + 1;


	    if (FindCommonElements(myTriangles[Tri_org[t].id[i]], myTriangles[Tri_org[t].id[j]], common)){
		if (common[0] != 2){
		    ErrWarnMessage(__LINE__, "CMeshImp::FindNeighbourTriangles:: Non-manifold surface", 0);
		}
		else{
		    Tri_org[t].neighbour[i] = (common[1] == t) ? common[2] : common[1];
		}
	    }
	    else{
		ErrWarnMessage(__LINE__, "CMeshImp::FindNeighbourTriangles::Non-manifold surface", 0);
	    }
	}
    }


    //clean up
    //for (int i = 0; i < numVert_org; i++){ delete[] myTriangles[i]; }

}
void CMeshImp::InitialSortAroundVertex(int**myTriangles)
{
    //Sort the triangle fan around each vertex
    //only used in the initilization stage since it depends of the input triangulation

    int common_tri[20];
    for (int ip = 0; ip < numVert_imp; ip++){
	int iq_prv;
	for (int i = 1; i < Vert_imp[ip].connect[0]; i++){
	    int iq = Vert_imp[ip].connect[i];
	    //get the two triangles that share the edge ip-iq
	    FindCommonElements(myTriangles[ip], myTriangles[iq],common_tri);
	    if (common_tri[0] != 2){
		ErrWarnMessage(__LINE__, "CMeshImp::InitialSortAroundVertex:: Input mesh is not watertight", 0);
	    }

	    int tri1(common_tri[1]), tri2(common_tri[2]), ik1, ik2;
	    //ik1 and ik2 are the third vertex in tr1 and tri2
	    for (int j = 0; j < 3; j++){
		if (Tri_org[tri1].id[j] != ip &&Tri_org[tri1].id[j] != iq){ ik1 = Tri_org[tri1].id[j]; }
	    }
	    for (int j = 0; j < 3; j++){
		if (Tri_org[tri2].id[j] != ip &&Tri_org[tri2].id[j] != iq){ ik2 = Tri_org[tri2].id[j]; }
	    }

	    int ik; //the actual replacement
	    if (i == 1){ ik = ik1; }//anyone of them would fit (this can be further utilized to make the sorting consistent i.e., CCW or CW)
	    else{
		//check
		if (iq_prv != ik1&&iq_prv != ik2){
		    ErrWarnMessage(__LINE__, "CMeshImp::InitialSortAroundVertex:: Input mesh has invalid connnectivity #1", 0);
		}
		ik = (iq_prv == ik1) ? ik2 : ik1;
	    }

	    int  ik_id = GetIndex(ik, Vert_imp[ip].connect);

	    if (ik_id < 0){
		ErrWarnMessage(__LINE__, "CMeshImp::InitialSortAroundVertex:: Input mesh has invalid connnectivity #2", 0);
	    }

	    std::swap(Vert_org[ip].connect[ik_id], Vert_org[ip].connect[i + 1]);
	    std::swap(Vert_imp[ip].connect[ik_id], Vert_imp[ip].connect[i + 1]);



	    iq_prv = iq;
	}

    }
}
void CMeshImp::ScaleInputMesh()
{
    //scale input mesh inside the unit box (0,0,0) (1,1,1)

    for (int i = 0; i < numVert_org; i++){
	Vert_org[i].x[0] = Vert_imp[i].x[0] = Vert_imp[i].x[0] - myBoundingBox.xmin;
	Vert_org[i].x[1] = Vert_imp[i].x[1] = Vert_imp[i].x[1] - myBoundingBox.ymin;
	Vert_org[i].x[2] = Vert_imp[i].x[2] = Vert_imp[i].x[2] - myBoundingBox.ymin;
    }

    scale_factor = 0.0;
    scale_factor = std::max(myBoundingBox.lx, myBoundingBox.ly);
    scale_factor = std::max(scale_factor, myBoundingBox.lz);

    for (int i = 0; i < numVert_org; i++){
	Vert_org[i].x[0] = Vert_imp[i].x[0] = Vert_imp[i].x[0] / scale_factor;
	Vert_org[i].x[1] = Vert_imp[i].x[1] = Vert_imp[i].x[1] / scale_factor;
	Vert_org[i].x[2] = Vert_imp[i].x[2] = Vert_imp[i].x[2] / scale_factor;
    }

    myBoundingBox.xmax -= myBoundingBox.xmin; myBoundingBox.ymax -= myBoundingBox.ymin; myBoundingBox.zmax -= myBoundingBox.zmin;
    myBoundingBox.xmax /= scale_factor; myBoundingBox.ymax /= scale_factor; myBoundingBox.zmax /= scale_factor;
    myBoundingBox.xmin = 0.0; myBoundingBox.ymin = 0.0; myBoundingBox.zmin = 0.0;
    myBoundingBox.lx = myBoundingBox.xmax - myBoundingBox.xmin; myBoundingBox.ly = myBoundingBox.ymax - myBoundingBox.ymin; myBoundingBox.lz = myBoundingBox.zmax - myBoundingBox.zmin;

}
void CMeshImp::GetFanTriangle()
{
    for (int t = 0; t < numTri_org; t++){
	for (int i = 0; i < 3; i++){
	    int j = (i == 2) ? 0 : i + 1;
	    int vert_i(Tri_org[t].id[i]), vert_j(Tri_org[t].id[j]);
	    int entry_index = GetIndex(vert_i, Vert_imp[vert_j].connect);
	    if (entry_index == -1){
		Vert_imp[vert_i].connect[++Vert_imp[vert_i].connect[0]] = Vert_org[vert_i].connect[++Vert_org[vert_i].connect[0]] = vert_j;
		Vert_imp[vert_j].connect[++Vert_imp[vert_j].connect[0]] = Vert_org[vert_j].connect[++Vert_org[vert_j].connect[0]] = vert_i;
	    }

	    if (Vert_imp[vert_i].connect[0] >= MAX_CONNECT || Vert_imp[vert_j].connect[0] >= MAX_CONNECT){
		ErrWarnMessage(__LINE__, "CMeshImp::CMeshImp:: increase MAX_CONNECT in Common.h", 0);
	    }
	}
    }
}


//** Sifting/Simplification **//
void CMeshImp::Simp(int targetNumSamples, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose)
{
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "****************************** Simplification Started ***************************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    Statistics myStats;
    myStats.GetAllStats(numVert_imp, Vert_imp);

    double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;
    Operator SimpOpt(Vert_org, &isObtuse, &isAcute);
    SimpOpt.constraints.isSmooth = isSmooth;
    SimpOpt.constraints.dev = devFactor;
    SimpOpt.constraints.isDelaunay = isDelaunay;

    SimpOpt.constraints.isMinAngle = true;
    SimpOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;

    SimpOpt.constraints.isMaxAngle = true;
    SimpOpt.constraints.MaxAngle = (maxAngleAllow<0.0) ? myStats.GetMaxAngle() : maxAngleAllow;

    SimpOpt.constraints.isNonobtuse = false;
    SimpOpt.constraints.isEdgeLen = false;
    SimpOpt.constraints.isMaximal = false;


    int verticesHandler[4];
    verticesHandler[0] = 2;


    int itter = 0;
    double vertex_void[3];

    if (verbose){
	fprintf(stdout, "\nRemoving tri-valent vertices --> ");
    }
    int numVert_before = numVert_imp;
    auto start_time = std::chrono::high_resolution_clock::now();
    SimpOpt.TriValentRemoval(numVert_imp, Vert_imp);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (verbose){
	fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);

	myStats.GetAllStats(numVert_imp, Vert_imp);
	DisplayStat(numVert_imp, Vert_imp,0);
	fprintf(stdout, "\n itter: %i",itter);
    }


    std::chrono::duration<double> elapsed_time_acc;


    while (true){
	numVert_before = numVert_imp;
	itter++;

	auto start_time1 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < numVert_imp; i++){
	    indices.push_back(i);
	}
	random_shuffle(indices.begin(), indices.end());


	//******* Ejection Two
	verticesHandler[0] = 2;
	//for (int ip1 = 0; ip1 < numVert_imp; ip1++)
	for (int id = 0; id < numVert_imp; id++){
	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp){ continue; }


	    if (Vert_imp[ip1].connect[0] == 0) { continue; }

	    bool ejected = false;

	    verticesHandler[1] = ip1;
	    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

		verticesHandler[2] = Vert_imp[ip1].connect[i];

		vertex_void[0] = (Vert_imp[verticesHandler[1]].x[0] +
			Vert_imp[verticesHandler[2]].x[0])/ 2.0;

		vertex_void[1] = (Vert_imp[verticesHandler[1]].x[1] +
			Vert_imp[verticesHandler[2]].x[1]) / 2.0;

		vertex_void[2] = (Vert_imp[verticesHandler[1]].x[2] +
			Vert_imp[verticesHandler[2]].x[2]) / 2.0;


		int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
		if (SimpOpt.constraints.dev > 0){
		    SimpOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - SimpOpt.constraints.dev;
		}

		if (SimpOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
		    ejected = true;
		    break;
		}
	    }

	    if (numVert_imp <= targetNumSamples){ break; }
	}

	if (numVert_imp > targetNumSamples){

	    //******* Ejection Three
	    verticesHandler[0] = 3;
	    for (int id = 0; id < numVert_imp; id++){
		int ip1 = indices[id];
		if (ip1 >= numVert_imp){ continue; }

		if (Vert_imp[ip1].connect[0] == 0) { continue; }

		verticesHandler[1] = ip1;

		verticesHandler[3] = Vert_imp[ip1].connect[Vert_imp[ip1].connect[0]];

		for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

		    verticesHandler[2] = Vert_imp[ip1].connect[i];

		    vertex_void[0] = (Vert_imp[verticesHandler[1]].x[0] +
			    Vert_imp[verticesHandler[2]].x[0] +
			    Vert_imp[verticesHandler[3]].x[0]) / 3.0;

		    vertex_void[1] = (Vert_imp[verticesHandler[1]].x[1] +
			    Vert_imp[verticesHandler[2]].x[1] +
			    Vert_imp[verticesHandler[3]].x[1]) / 3.0;

		    vertex_void[2] = (Vert_imp[verticesHandler[1]].x[2] +
			    Vert_imp[verticesHandler[2]].x[2] +
			    Vert_imp[verticesHandler[3]].x[2]) / 3.0;


		    int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
		    if (SimpOpt.constraints.dev > 0){
			SimpOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - SimpOpt.constraints.dev;
		    }

		    if (SimpOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			break;
		    }

		    verticesHandler[3] = verticesHandler[2];
		}

		if (numVert_imp <= targetNumSamples){ break; }
	    }
	}

	auto end_time1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

	elapsed_time_acc = elapsed_time_acc + elapsed_time1;


	if (verbose){
	    myStats.GetAllStats(numVert_imp, Vert_imp);

	    DisplayStat(numVert_imp, Vert_imp,0);

	    fprintf(stdout, "\n Reducation Ratio= %f %", 100.0*(double(numVert_org - numVert_imp) / double(numVert_org)));

	    fprintf(stdout, "\n itter: %i", itter);
	}

	if (numVert_imp <= targetNumSamples){ break; }
	if (numVert_before == numVert_imp){ break; }

    }

    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "***************************** Simplification Finished ***************************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    std::cout << " \nTotal elapsed time: " << elapsed_time_acc.count() << " (s)\n";

    GetMeshOBJ("output.obj", 0, "");
}

//** Non-obtuse Remeshing **//
void CMeshImp::NonobtuseRemeshing(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose)
{
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "************************ Non-obtuse Remeshing Started ***************************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;
    isObtuse = false;//no partial improvement
    isAcute = false;

    Statistics myStats;
    myStats.GetAllStats(numVert_imp, Vert_imp);


    Operator NonObtuseOpt(Vert_org, &isObtuse, &isAcute);

    NonObtuseOpt.constraints.isSmooth = isSmooth;
    NonObtuseOpt.constraints.dev = devFactor;
    NonObtuseOpt.constraints.isDelaunay = isDelaunay;

    NonObtuseOpt.constraints.isMinAngle = true;
    NonObtuseOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;
    NonObtuseOpt.constraints.isMaxAngle = true;
    NonObtuseOpt.constraints.MaxAngle = 90.0;

    NonObtuseOpt.constraints.isNonobtuse = true;

    NonObtuseOpt.constraints.isEdgeLen = false;
    NonObtuseOpt.constraints.isMaximal = false;


    int verticesHandler[4];
    verticesHandler[0] = 2;

    int ip2, ip3;
    int itter = 0;
    double vertex_void[3];

    if (verbose){
	fprintf(stdout, "\nRemoving tri-valent vertices --> ");
    }

    int numVert_before = numVert_imp;

    auto start_time = std::chrono::high_resolution_clock::now();
    NonObtuseOpt.TriValentRemoval(numVert_imp, Vert_imp);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (verbose){
	fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
    }

    myStats.GetAllStats(numVert_imp, Vert_imp);

    if (verbose){
	DisplayStat(numVert_imp, Vert_imp, 1);
	fprintf(stdout, "\n itter: %i", itter);
    }


    std::chrono::duration<double> elapsed_time_acc = end_time-start_time;

    while (myStats.GetNumNonObtuse() > 0){
	itter++;

	auto start_time1 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < numVert_imp; i++){
	    indices.push_back(i);
	}
	random_shuffle(indices.begin(), indices.end());

	//******* Relocation
	for (size_t id = 0; id < indices.size(); id++){
	    int i = indices[id];
	    if (i >= numVert_imp){ continue; }

	    //printf("\nRelocating [%d]", i);
	    //fflush(stdout);


	    int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[i].x);

	    //only do smoothness at smooth areas
	    //i.e., areas that deviates from perfect smooth (180.0 deg) by amount equal to dev factor
	    //the area is defined by closest surface vertex to i

	    if (NonObtuseOpt.constraints.dev > 0){
		NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
	    }

	    //if (ObtuseHead(i, ip2, ip3)){
	    NonObtuseOpt.Relocation(i, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);
	    //}

	}


	//******* Ejection
	verticesHandler[0] = 2;
	for (size_t id = 0; id < indices.size(); id++){
	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp){ continue; }


	    if (Vert_imp[ip1].connect[0] == 0) { continue; }
	    if (ObtuseHead(ip1, ip2, ip3)){
		//we could eject either ip1, ip2 or ip3 (along with another vertex)
		bool ejected = false;

		verticesHandler[1] = ip1;
		for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

		    verticesHandler[2] = Vert_imp[ip1].connect[i];
		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;
		    int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
		    if (NonObtuseOpt.constraints.dev >0){
			NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
		    }

		    if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			ejected = true;
			break;
		    }
		}

		if (!ejected){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
			int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
			if (NonObtuseOpt.constraints.dev >0){
			    NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
			}

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    ejected = true;
			    break;
			}
		    }
		}

		if (!ejected){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
			int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
			if (NonObtuseOpt.constraints.dev >0){
			    NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
			}


			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    ejected = true;
			    break;
			}
		    }
		}
	    }
	}


	//******* Injection
	verticesHandler[0] = 3;
	for (size_t id = 0; id < indices.size(); id++){
	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp){ continue; }

	    if (Vert_imp[ip1].connect[0] == 0) { continue; }
	    if (ObtuseHead(ip1, ip2, ip3)){

		//call operator on triangle (ip1, ip2,ip3)
		verticesHandler[1] = ip1;
		verticesHandler[2] = ip2;
		verticesHandler[3] = ip3;

		vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		int closestSurfID = surfaceKdTree.FindNearest(vertex_void);

		if (NonObtuseOpt.constraints.dev >0){
		    NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
		}


		NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0);
	    }
	}

	//******* Attractor Ejection
	verticesHandler[0] = 2;
	for (size_t id = 0; id < indices.size(); id++){
	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp){ continue; }

	    if (Vert_imp[ip1].connect[0] == 0) { continue; }
	    if (ObtuseHead(ip1, ip2, ip3)){
		//we could eject either ip1, ip2 or ip3 (along with another vertex)
		bool ejected = false;

		verticesHandler[1] = ip1;
		for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){

		    verticesHandler[2] = Vert_imp[ip1].connect[i];
		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
		    if (NonObtuseOpt.constraints.dev >0){
			NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
		    }

		    if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			ejected = true;
			break;
		    }
		}

		if (!ejected){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

			int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
			if (NonObtuseOpt.constraints.dev >0){
			    NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
			}

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    ejected = true;
			    break;
			}
		    }
		}

		if (!ejected){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
			int closestSurfID = surfaceKdTree.FindNearest(vertex_void);
			if (NonObtuseOpt.constraints.dev >0){
			    NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
			}

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    ejected = true;
			    break;
			}
		    }
		}
	    }
	}

	//******* Repeller Injection
	verticesHandler[0] = 3;
	for (size_t id = 0; id < indices.size(); id++){
	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp){ continue; }

	    if (Vert_imp[ip1].connect[0] == 0) { continue; }
	    if (ObtuseHead(ip1, ip2, ip3)){

		//call operator on triangle (ip1, ip2,ip3)
		verticesHandler[1] = ip1;
		verticesHandler[2] = ip2;
		verticesHandler[3] = ip3;

		vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		int closestSurfID = surfaceKdTree.FindNearest(vertex_void);

		if (NonObtuseOpt.constraints.dev >0){
		    NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
		}

		NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);

	    }
	}

	auto end_time1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

	elapsed_time_acc = elapsed_time_acc + elapsed_time1;

	myStats.GetAllStats(numVert_imp, Vert_imp);

	if (verbose){
	    DisplayStat(numVert_imp, Vert_imp,1);
	    fprintf(stdout, "\n itter: %i", itter);
	}
    }

    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "************************ Non-obtuse Remeshing Finished **************************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";

    GetMeshOBJ("output.obj", 0, "");
}

void CMeshImp::NonobtuseRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, bool verbose)
{
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "******************* Non-obtuse Remeshing - Interleave Started *******************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;

    Statistics myStats;
    myStats.GetAllStats(numVert_imp, Vert_imp);

    //Set the constriants handler
    Operator NonObtuseOpt(Vert_org, &isObtuse, &isAcute);

    NonObtuseOpt.constraints.isSmooth = isSmooth;
    NonObtuseOpt.constraints.dev = devFactor;
    NonObtuseOpt.constraints.isDelaunay = isDelaunay;

    NonObtuseOpt.constraints.isMinAngle = true;
    NonObtuseOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;
    NonObtuseOpt.constraints.isMaxAngle = true;
    NonObtuseOpt.constraints.MaxAngle = 90;

    NonObtuseOpt.constraints.isNonobtuse = true;

    NonObtuseOpt.constraints.isEdgeLen = false;
    NonObtuseOpt.constraints.isMaximal = false;


    int verticesHandler[4];


    int ip2, ip3;
    int itter = 0;
    double vertex_void[3];

    if (verbose){
	fprintf(stdout, "\nRemoving tri-valent vertices --> ");
    }
    int numVert_before = numVert_imp;

    auto start_time = std::chrono::high_resolution_clock::now();
    NonObtuseOpt.TriValentRemoval(numVert_imp, Vert_imp);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (verbose){
	fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
    }


    myStats.GetAllStats(numVert_imp, Vert_imp);

    if (verbose){
	DisplayStat(numVert_imp, Vert_imp,1);
	fprintf(stdout, "\n itter: %i", itter);
    }

    std::chrono::duration<double> elapsed_time_acc = end_time - start_time;

    while (myStats.GetNumNonObtuse() > 0){
	itter++;

	auto start_time1 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < numVert_imp; i++){
	    indices.push_back(i);
	}
	random_shuffle(indices.begin(), indices.end());

	// Relocation
	for (size_t id = 0; id < indices.size(); id++){

	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp || Vert_imp[ip1].connect[0] == 0){ continue; }


	    int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[ip1].x);

	    if (NonObtuseOpt.constraints.dev >0){
		NonObtuseOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - NonObtuseOpt.constraints.dev;
	    }

	    isObtuse = ObtuseHead(ip1, ip2, ip3);
	    //Relocation
	    NonObtuseOpt.Relocation(ip1, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);


	    //test if it is still an obtuse head
	    if (ObtuseHead(ip1, ip2, ip3)){
		isObtuse = true;

		bool fixed = false;

		//Ejection
		verticesHandler[0] = 2;

		if (!fixed){
		    verticesHandler[1] = ip1;
		    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip1].connect[i];
			vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    //Injection
		    verticesHandler[0] = 3;
		    verticesHandler[1] = ip1;
		    verticesHandler[2] = ip2;
		    verticesHandler[3] = ip3;

		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    if (NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0)){
			fixed = true;
		    }
		}

		//Attracotr Ejection
		verticesHandler[0] = 2;
		if (!fixed){
		    verticesHandler[1] = ip1;
		    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip1].connect[i];
			vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
			vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
			vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}
		if (!fixed){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

			if (NonObtuseOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}

		//Repeller Injection
		if (!fixed){
		    verticesHandler[0] = 3;
		    verticesHandler[1] = ip1;
		    verticesHandler[2] = ip2;
		    verticesHandler[3] = ip3;

		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    NonObtuseOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);
		}
	    }
	}

	auto end_time1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

	elapsed_time_acc = elapsed_time_acc + elapsed_time1;

	myStats.GetAllStats(numVert_imp, Vert_imp);

	if (verbose){
	    DisplayStat(numVert_imp, Vert_imp,1);
	    fprintf(stdout, "\n itter: %i", itter);
	}
    }

    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "****************** Non-obtuse Remeshing - Interleave Finished *******************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";

    GetMeshOBJ("output.obj", 0, "");
}
void CMeshImp::AcuteRemeshing_InterleaveOpt(int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double minAngleAllow, double maxAngleAllow, bool verbose)
{
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "********************** Acute Remeshing - Interleave Started *********************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;

    Statistics myStats;
    myStats.GetAllStats(numVert_imp, Vert_imp);

    //Set the constriants handler
    Operator AcuteOpt(Vert_org, &isObtuse, &isAcute);

    AcuteOpt.constraints.isSmooth = isSmooth;
    AcuteOpt.constraints.dev = devFactor;
    AcuteOpt.constraints.isDelaunay = isDelaunay;

    AcuteOpt.constraints.isMinAngle = true;
    AcuteOpt.constraints.MinAngle = (minAngleAllow<0.0) ? myStats.GetMinAngle() : minAngleAllow;
    AcuteOpt.constraints.isMaxAngle = true;
    AcuteOpt.constraints.MaxAngle = maxAngleAllow;

    AcuteOpt.constraints.isNonobtuse = false;

    AcuteOpt.constraints.isEdgeLen = false;
    AcuteOpt.constraints.isMaximal = false;


    int verticesHandler[4];


    int ip2, ip3;
    int itter = 0;
    double vertex_void[3];

    if (verbose){
	fprintf(stdout, "\nRemoving tri-valent vertices --> ");
    }

    int numVert_before = numVert_imp;

    auto start_time = std::chrono::high_resolution_clock::now();
    AcuteOpt.TriValentRemoval(numVert_imp, Vert_imp);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (verbose){
	fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
    }


    myStats.GetAllStats(numVert_imp, Vert_imp);
    int numAcute = myStats.CalcNumAcute(numVert_imp, Vert_imp, maxAngleAllow);

    if (verbose){
	DisplayStat(numVert_imp, Vert_imp,1);
	fprintf(stdout, "\n #Acute triangle = %i", numAcute);
	fprintf(stdout, "\n itter: %i", itter);
    }


    std::chrono::duration<double> elapsed_time_acc = end_time - start_time;

    while (numAcute > 0){
	itter++;

	auto start_time1 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < numVert_imp; i++){
	    indices.push_back(i);
	}
	random_shuffle(indices.begin(), indices.end());

	// Relocation
	for (size_t id = 0; id < indices.size(); id++){
	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp || Vert_imp[ip1].connect[0] == 0){ continue; }


	    int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[ip1].x);

	    if (AcuteOpt.constraints.dev > 0){
		AcuteOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - AcuteOpt.constraints.dev;
	    }

	    isAcute = AcuteHead(ip1, ip2, ip3, maxAngleAllow);

	    //Relocation
	    AcuteOpt.Relocation(ip1, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);

	    //test if it is still an obtuse head
	    if (AcuteHead(ip1, ip2, ip3, maxAngleAllow)){
		isAcute = false;


		bool fixed = false;

		//Ejection
		verticesHandler[0] = 2;

		if (!fixed){
		    verticesHandler[1] = ip1;
		    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip1].connect[i];
			vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;
			if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
			if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
			if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    //Injection
		    verticesHandler[0] = 3;
		    verticesHandler[1] = ip1;
		    verticesHandler[2] = ip2;
		    verticesHandler[3] = ip3;

		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    if (AcuteOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0)){
			fixed = true;
		    }
		}

		//Attracotr Ejection
		verticesHandler[0] = 2;
		if (!fixed){
		    verticesHandler[1] = ip1;
		    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip1].connect[i];
			vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
			vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
			vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;
			if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;
			if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}
		if (!fixed){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;
			if (AcuteOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}

		//Repeller Injection
		if (!fixed){
		    verticesHandler[0] = 3;
		    verticesHandler[1] = ip1;
		    verticesHandler[2] = ip2;
		    verticesHandler[3] = ip3;

		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    AcuteOpt.Injection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);
		}
	    }
	}

	auto end_time1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

	elapsed_time_acc = elapsed_time_acc + elapsed_time1;

	myStats.GetAllStats(numVert_imp, Vert_imp);
	numAcute = myStats.CalcNumAcute(numVert_imp, Vert_imp, maxAngleAllow);

	if (verbose){
	    DisplayStat(numVert_imp, Vert_imp, 1);
	    fprintf(stdout, "\n #Acute triangle = %i", numAcute);
	    fprintf(stdout, "\n itter: %i", itter);
	}
    }

    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "********************** Acute Remeshing - Interleave Finished ********************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";

    GetMeshOBJ("output.obj", 0, "");
}

//** Small Angle Elimination **//
void CMeshImp::SmallAngleElimination_InterleaveOpt(double targetMinAngle, int samplingBudget, int numSurfaceLayer, bool isSmooth, double theta_d, bool isDelaunay, double maxAngleAllow, bool verbose)
{
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "****************** Small Angle Elimination - Interleave Started *****************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    double devFactor = (theta_d < 0.0) ? -1.0 : 180 - theta_d;

    Statistics myStats;
    myStats.GetAllStats(numVert_imp, Vert_imp);

    //Set the constriants handler
    Operator SmallAngOpt (Vert_org, &isObtuse, &isAcute);

    SmallAngOpt.constraints.isSmooth = isSmooth;
    SmallAngOpt.constraints.dev = devFactor;
    SmallAngOpt.constraints.isDelaunay = isDelaunay;

    SmallAngOpt.constraints.isMinAngle = true;
    SmallAngOpt.constraints.MinAngle = std::max(2.0, 1.5*myStats.GetMinAngle());

    SmallAngOpt.constraints.isMaxAngle = true;
    SmallAngOpt.constraints.MaxAngle = (maxAngleAllow<0.0) ? myStats.GetMaxAngle() : maxAngleAllow;

    SmallAngOpt.constraints.isNonobtuse = false;

    SmallAngOpt.constraints.isEdgeLen = false;
    SmallAngOpt.constraints.isMaximal = false;


    int verticesHandler[4];


    int ip2, ip3;
    int itter = 0;
    double vertex_void[3];

    if (verbose){
	fprintf(stdout, "\nRemoving tri-valent vertices --> ");
    }

    int numVert_before = numVert_imp;

    auto start_time = std::chrono::high_resolution_clock::now();
    SmallAngOpt.TriValentRemoval(numVert_imp, Vert_imp);
    auto end_time = std::chrono::high_resolution_clock::now();

    if (verbose){
	fprintf(stdout, "%d vertex removed\n", numVert_before - numVert_imp);
    }

    myStats.GetAllStats(numVert_imp, Vert_imp);

    if (verbose){
	DisplayStat(numVert_imp, Vert_imp,0);
	fprintf(stdout, "\n itter: %i", itter);
    }

    std::chrono::duration<double> elapsed_time_acc = end_time - start_time;
    isObtuse = false;

    while (true){
	itter++;

	auto start_time1 = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < numVert_imp; i++){
	    indices.push_back(i);
	}
	random_shuffle(indices.begin(), indices.end());

	// Relocation
	for (size_t id = 0; id < indices.size(); id++){

	    int ip1 = indices[id];
	    if (ip1 >= numVert_imp || Vert_imp[ip1].connect[0] == 0){ continue; }


	    //test if it is still an obtuse head
	    if (TooSmallHeadAngle(ip1, ip2, ip3, targetMinAngle)){

		int closestSurfID = surfaceKdTree.FindNearest(Vert_imp[ip1].x);

		if (SmallAngOpt.constraints.dev >0){
		    SmallAngOpt.constraints.isSmooth = Vert_org[closestSurfID].dih_ang > 180.0 - SmallAngOpt.constraints.dev;
		}

		bool fixed = false;

		//Relocation
		fixed = SmallAngOpt.Relocation(ip1, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer);


		//Ejection
		verticesHandler[0] = 2;

		if (!fixed){
		    verticesHandler[1] = ip1;
		    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip1].connect[i];
			vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[Vert_imp[ip1].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[Vert_imp[ip1].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[Vert_imp[ip1].connect[i]].x[2]) / 2.0;


			if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

			if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

			if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 0)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    //Injection
		    verticesHandler[0] = 3;
		    verticesHandler[1] = ip1;
		    verticesHandler[2] = ip2;
		    verticesHandler[3] = ip3;

		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    if (SmallAngOpt.AggressiveInjection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 0)){
			fixed = true;
		    }
		}

		//Attracotr Ejection
		verticesHandler[0] = 2;
		if (!fixed){
		    verticesHandler[1] = ip1;
		    for (int i = 1; i <= Vert_imp[ip1].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip1].connect[i];
			vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
			vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
			vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

			if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}

		if (!fixed){
		    verticesHandler[1] = ip2;
		    for (int i = 1; i <= Vert_imp[ip2].connect[0]; i++){

			verticesHandler[2] = Vert_imp[ip2].connect[i];
			vertex_void[0] = (Vert_imp[ip2].x[0] + Vert_imp[Vert_imp[ip2].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip2].x[1] + Vert_imp[Vert_imp[ip2].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip2].x[2] + Vert_imp[Vert_imp[ip2].connect[i]].x[2]) / 2.0;

			if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}
		if (!fixed){
		    verticesHandler[1] = ip3;
		    for (int i = 1; i <= Vert_imp[ip3].connect[0]; i++){
			verticesHandler[2] = Vert_imp[ip3].connect[i];
			vertex_void[0] = (Vert_imp[ip3].x[0] + Vert_imp[Vert_imp[ip3].connect[i]].x[0]) / 2.0;
			vertex_void[1] = (Vert_imp[ip3].x[1] + Vert_imp[Vert_imp[ip3].connect[i]].x[1]) / 2.0;
			vertex_void[2] = (Vert_imp[ip3].x[2] + Vert_imp[Vert_imp[ip3].connect[i]].x[2]) / 2.0;

			if (SmallAngOpt.Ejection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, samplingBudget, numSurfaceLayer, 1)){
			    fixed = true;
			    break;
			}
		    }
		}

		//Repeller Injection
		if (!fixed){
		    verticesHandler[0] = 3;
		    verticesHandler[1] = ip1;
		    verticesHandler[2] = ip2;
		    verticesHandler[3] = ip3;

		    vertex_void[0] = (Vert_imp[ip1].x[0] + Vert_imp[ip2].x[0] + Vert_imp[ip3].x[0]) / 3.0;
		    vertex_void[1] = (Vert_imp[ip1].x[1] + Vert_imp[ip2].x[1] + Vert_imp[ip3].x[1]) / 3.0;
		    vertex_void[2] = (Vert_imp[ip1].x[2] + Vert_imp[ip2].x[2] + Vert_imp[ip3].x[2]) / 3.0;

		    SmallAngOpt.AggressiveInjection(verticesHandler, numVert_imp, Vert_imp, closestSurfID, 1, numSurfaceLayer, 1);
		}
	    }
	}

	auto end_time1 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_time1 = end_time1 - start_time1;

	elapsed_time_acc = elapsed_time_acc + elapsed_time1;

	myStats.GetAllStats(numVert_imp, Vert_imp);

	if (verbose){
	    DisplayStat(numVert_imp, Vert_imp,0);
	    fprintf(stdout, "\n itter: %i", itter);
	}

	if (myStats.GetMinAngle() > targetMinAngle){ break; }


	SmallAngOpt.constraints.MinAngle = std::max(2.0, 1.5*myStats.GetMinAngle());
    }

    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
    fprintf(stdout, "****************** Small Angle Elimination - Interleave Finished ****************");
    fprintf(stdout, "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");

    std::cout << "\n Total elapsed time: " << elapsed_time_acc.count() << " (s)\n";

    GetMeshOBJ("output.obj", 0, "");
}

//** Evaluate head **//
bool CMeshImp::ObtuseHead(int ip, int&ip1, int&ip2)
{

    ip2 = Vert_imp[ip].connect[Vert_imp[ip].connect[0]];

    for (int i = 1; i <= Vert_imp[ip].connect[0]; i++){
	ip1 = Vert_imp[ip].connect[i];

	double angle = AngleVectVect(Vert_imp[ip1].x[0] - Vert_imp[ip].x[0], Vert_imp[ip1].x[1] - Vert_imp[ip].x[1], Vert_imp[ip1].x[2] - Vert_imp[ip].x[2],
		Vert_imp[ip2].x[0] - Vert_imp[ip].x[0], Vert_imp[ip2].x[1] - Vert_imp[ip].x[1], Vert_imp[ip2].x[2] - Vert_imp[ip].x[2])*RadToDeg; /* 57.295779513078550 = 180.0 / PI*/

	if (angle > 90.0 + _tol ){
	    return true;
	}
	ip2 = ip1;
    }

    return false;

}
bool CMeshImp::AcuteHead(int ip, int&ip1, int&ip2, double measureAngle)
{

    ip2 = Vert_imp[ip].connect[Vert_imp[ip].connect[0]];

    for (int i = 1; i <= Vert_imp[ip].connect[0]; i++){
	ip1 = Vert_imp[ip].connect[i];

	double angle = AngleVectVect(Vert_imp[ip1].x[0] - Vert_imp[ip].x[0], Vert_imp[ip1].x[1] - Vert_imp[ip].x[1], Vert_imp[ip1].x[2] - Vert_imp[ip].x[2],
		Vert_imp[ip2].x[0] - Vert_imp[ip].x[0], Vert_imp[ip2].x[1] - Vert_imp[ip].x[1], Vert_imp[ip2].x[2] - Vert_imp[ip].x[2])*RadToDeg;

	if (angle > measureAngle + _tol){
	    return true;
	}
	ip2 = ip1;
    }

    return false;

}
bool CMeshImp::TooSmallHeadAngle(int ip, int&ip1, int&ip2, double measureAngle)
{

    ip2 = Vert_imp[ip].connect[Vert_imp[ip].connect[0]];

    for (int i = 1; i <= Vert_imp[ip].connect[0]; i++){
	ip1 = Vert_imp[ip].connect[i];

	double angle = AngleVectVect(Vert_imp[ip1].x[0] - Vert_imp[ip].x[0], Vert_imp[ip1].x[1] - Vert_imp[ip].x[1], Vert_imp[ip1].x[2] - Vert_imp[ip].x[2],
		Vert_imp[ip2].x[0] - Vert_imp[ip].x[0], Vert_imp[ip2].x[1] - Vert_imp[ip].x[1], Vert_imp[ip2].x[2] - Vert_imp[ip].x[2])*RadToDeg;

	if (angle + _tol< measureAngle){
	    return true;
	}
	ip2 = ip1;
    }

    return false;

}

//** Postprocessing **//
void CMeshImp::WriteStatToFile(std::string filename, int numV, vert*Vert, bool obt)
{
    Statistics myStats;
    myStats.GetAllStats(numV, Vert);

    FILE *fp = fopen(filename.c_str(), "w");
    myStats.DisplayStatistics(fp, obt);
    fclose(fp);

}
void CMeshImp::DisplayStat(int numV, vert*Vert, bool obt)
{
    Statistics myStats;
    myStats.GetAllStats(numV, Vert);
    myStats.DisplayStatistics(stdout, obt);

}
void CMeshImp::GetMeshOBJ(std::string filename, bool obtuse, std::string filename_obt)
{
    //write mesh in .obj file
    fprintf(stdout, "\n Writing output file.\n");

    std::fstream file(filename.c_str(), std::ios::out);
    std::fstream file_obt(filename_obt.c_str(), std::ios::out);



    file << "#V " <<numVert_imp << std::endl;
    for (int id0 = 0; id0 < numVert_imp; id0++){
	file << "v " << Vert_imp[id0].x[0] << " " << Vert_imp[id0].x[1] << " " << Vert_imp[id0].x[2] << std::endl;
	if (obtuse){
	    file_obt << "v " << Vert_imp[id0].x[0] << " " << Vert_imp[id0].x[1] << " " << Vert_imp[id0].x[2] << std::endl;
	}
    }

    for (int id0 = 0; id0 < numVert_imp; id0++){
	int id1 = Vert_imp[id0].connect[Vert_imp[id0].connect[0]];

	for (int i = 1; i <= Vert_imp[id0].connect[0]; i++){
	    int id2 = Vert_imp[id0].connect[i];

	    if (id0 < id1 && id0 < id2){
		file << "f " << id0 + 1 << "  " << id1 + 1 << "  " << id2 + 1 << std::endl;
		if (obtuse){
		    double angle1 = AngleVectVect(Vert_imp[id1].x[0] - Vert_imp[id0].x[0], Vert_imp[id1].x[1] - Vert_imp[id0].x[1], Vert_imp[id1].x[2] - Vert_imp[id0].x[2],
			    Vert_imp[id2].x[0] - Vert_imp[id0].x[0], Vert_imp[id2].x[1] - Vert_imp[id0].x[1], Vert_imp[id2].x[2] - Vert_imp[id0].x[2])*RadToDeg;

		    double angle2 = AngleVectVect(Vert_imp[id0].x[0] - Vert_imp[id1].x[0], Vert_imp[id0].x[1] - Vert_imp[id1].x[1], Vert_imp[id0].x[2] - Vert_imp[id1].x[2],
			    Vert_imp[id2].x[0] - Vert_imp[id1].x[0], Vert_imp[id2].x[1] - Vert_imp[id1].x[1], Vert_imp[id2].x[2] - Vert_imp[id1].x[2])*RadToDeg;

		    double angle3 = AngleVectVect(Vert_imp[id1].x[0] - Vert_imp[id2].x[0], Vert_imp[id1].x[1] - Vert_imp[id2].x[1], Vert_imp[id1].x[2] - Vert_imp[id2].x[2],
			    Vert_imp[id0].x[0] - Vert_imp[id2].x[0], Vert_imp[id0].x[1] - Vert_imp[id2].x[1], Vert_imp[id0].x[2] - Vert_imp[id2].x[2])*RadToDeg;
		    if (angle2>90.0 + _tol || angle3>90.0 + _tol){
			file_obt << "f " << id0 + 1 << "  " << id1 + 1 << "  " << id2 + 1 << std::endl;
		    }
		}
	    }
	    id1 = id2;
	}
    }

    file.close();
}

void PointOnSphere(double&x1, double&y1, double&z1, double x_s, double y_s, double z_s, double r_s) // Setting point(x1,y1,z1) to be on the sphere of radius r_s centered at (x_s, y_s, z_s)
{
    double x_v(x1 - x_s), y_v(y1 - y_s), z_v(z1 - z_s);

    double n = sqrt(x_v*x_v + y_v*y_v + z_v*z_v);

    x1 = x_v*r_s / n + x_s;
    y1 = y_v*r_s / n + y_s;
    z1 = z_v*r_s / n + z_s;


}
void CMeshImp::PerturbVertices()
{
    //input is vertices one a sphere
    //move verices that are not head of obtuse angles
    //in random direction 10% of shorest edge its connected to

    //srand(time(NULL));
    int n1, n2;
    RndNum myRandNum;
    myRandNum.InitiateRndNum(rand());

    for (int ip = 0; ip < numVert_imp; ip++){
	if (!ObtuseHead(ip, n1, n2)){
	    double shortest_edge_len = 10E10;
	    for (int j = 1; j <= Vert_imp[ip].connect[0]; j++){
		int iq = Vert_imp[ip].connect[j];
		shortest_edge_len = std::min(shortest_edge_len, Dist(Vert_imp[ip].x[0], Vert_imp[ip].x[1], Vert_imp[ip].x[2],
			    Vert_imp[iq].x[0], Vert_imp[iq].x[1], Vert_imp[iq].x[2]));
	    }
	    shortest_edge_len = sqrt(shortest_edge_len);
	    shortest_edge_len /= 4.0;

	    Vert_imp[ip].x[0] += (myRandNum.RandNumGenerator())* shortest_edge_len;
	    Vert_imp[ip].x[1] += (myRandNum.RandNumGenerator())* shortest_edge_len;
	    Vert_imp[ip].x[2] += (myRandNum.RandNumGenerator())* shortest_edge_len;

	    PointOnSphere(Vert_imp[ip].x[0], Vert_imp[ip].x[1], Vert_imp[ip].x[2], 0.5, 0.5, 0.5, 0.5*sqrt(3.0));

	}
    }

    GetMeshOBJ("sphere8Perturbed.obj", 1, "input_obtuse.obj");
}

#endif //CMeshImpl

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8

