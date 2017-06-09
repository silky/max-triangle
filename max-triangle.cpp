#include "max-triangle.hpp"

namespace InscribedTriangle {
    bool is_convex(std::vector<mpq_class> polygon) {
	//check if polygon is a strictly convex polygon with vertices listed in CCW order

	unsigned int i;
	unsigned int sz = polygon.size()/2;
	mpq_class x;

	for (i = 0; i < sz; i++) {
	    x  = (polygon[(2*i+2)%(2*sz)] - polygon[(2*i+0)%(2*sz)])*(polygon[(2*i+5)%(2*sz)] - polygon[(2*i+1)%(2*sz)]);
	    x -= (polygon[(2*i+4)%(2*sz)] - polygon[(2*i+0)%(2*sz)])*(polygon[(2*i+3)%(2*sz)] - polygon[(2*i+1)%(2*sz)]);
	    if (x <= 0) return false;
	}

	return true;

    }


    void anchored_triangle(std::vector<mpq_class> polygon, mpq_class nx, mpq_class ny, unsigned int *reti, mpq_class *rett, state_flag *status){
	//find a triangle inscribed in a convex polygon that is of maximal area among all the ones with (nx,ny) as an outer normal
	//input:
	//  polygon: a list of coordinates (x_0, y_0, x_1, y_1, ..., x_{n-1}, y_{n-1}) of the vertices of a convex polygon
	//           the vertices must be in couter-clockwise order and all must be extreme points
	//  nx, ny : the desired normal vector
	//output:
	//  reti   : an array of 3 indices (i_a, i_b, i_c) between 0 and n-1
	//  rett   : an array of 2 rational numbers (s_b, s_c) between 0 and 1
	//           together, reti and rett specify the output triangle,
	//           z_a = z_{i_a}, z_b = z_{i_b}(1-s_b) + z_{i_b+1}s_b, and z_b = z_{i_b}(1-s_b) + z_{i_b+1}s_b
	//  status : a flag indicating the successful termination of the function or an error encountered

	mpq_class minx,maxx,x;
	unsigned int ai,bi,ci;
	mpq_class ax,ay,bx,by,cx,cy;
	mpq_class ebx,eby,ecx,ecy;
	mpq_class dbt, dct, t0, tq;
	mpq_class ee, neb, nec;
	mpq_class bt,ct;
	unsigned int i;
	unsigned int sz = polygon.size()/2;

	if (sz < 3) { *status = status_no_interior; return; } //a convex polygon has at least 3 vertices
	if (!is_convex(polygon)) { *status = status_not_convex; return; } //the vertices are not in convex position or not in counter-clockwise order
	*status = status_ok;

	maxx = minx = nx*polygon[0]+ny*polygon[1];
	ai = ci = 0;
	for (i = 1; i < sz; i++) {
	    x = nx*polygon[2*i+0]+ny*polygon[2*i+1];
	    if (x>maxx) { maxx=x; ci=i; }
	    if (x<minx) { minx=x; ai=i; }
	}
	ci+=sz;
	bi=ci-1;
	bt=1;
	ct=0;

	ax = polygon[(2*ai+0)%(2*sz)];
	ay = polygon[(2*ai+1)%(2*sz)];
	cx = polygon[(2*ci+0)%(2*sz)];
	cy = polygon[(2*ci+1)%(2*sz)];
	bx = cx;
	by = cy;
	ebx = polygon[(2*bi+2)%(2*sz)] - polygon[(2*bi+0)%(2*sz)];
	eby = polygon[(2*bi+3)%(2*sz)] - polygon[(2*bi+1)%(2*sz)];
	ecx = polygon[(2*ci+2)%(2*sz)] - polygon[(2*ci+0)%(2*sz)];
	ecy = polygon[(2*ci+3)%(2*sz)] - polygon[(2*ci+1)%(2*sz)];

	while(true) {
	    neb = ebx*nx + eby*ny; 
	    nec = ecx*nx + ecy*ny;
	    ee = ebx*ecy - ecx*eby;
	    if (neb < 0) { *status = status_runtime_error; return; }
	    if (nec > 0) { *status = status_runtime_error; return; }

	    if (neb == 0) { //eb.n = 0, advance b
		bi -= 1; 
		bt = 1;
		bx = polygon[(2*bi+2)%(2*sz)];
		by = polygon[(2*bi+3)%(2*sz)];
		ebx = polygon[(2*bi+2)%(2*sz)] - polygon[(2*bi+0)%(2*sz)];
		eby = polygon[(2*bi+3)%(2*sz)] - polygon[(2*bi+1)%(2*sz)];
		continue;
	    }
	    if (nec == 0) { //ec.n = 0, advance c
		ci += 1; 
		ct = 0;
		cx = polygon[(2*ci+0)%(2*sz)];
		cy = polygon[(2*ci+1)%(2*sz)];
		ecx = polygon[(2*ci+2)%(2*sz)] - polygon[(2*ci+0)%(2*sz)];
		ecy = polygon[(2*ci+3)%(2*sz)] - polygon[(2*ci+1)%(2*sz)];
		continue;
	    }
	    if (ee <= 0) { //eb and ec are parallel or converging, any advancement reduces area, so done
		reti[0] = ai; reti[1] = bi; reti[2] = ci;
		rett[0] = bt; rett[1] = ct;
		return;
	    }

	    tq  = nx*( (by-cy)*ebx*ecx + (cx-ax)*eby*ecx + (ax-bx)*ebx*ecy);
	    tq -= ny*( (bx-cx)*eby*ecy + (cy-ay)*ebx*ecy + (ay-by)*eby*ecx);
	    tq /= 2*ee;

	    dbt = -tq/neb;
	    dct = tq/nec;

	    if (dbt <= 0 || dct <= 0) {
		reti[0] = ai; reti[1] = bi; reti[2] = ci;
		rett[0] = bt; rett[1] = ct;
		return;
	    }

	    if (dbt < bt && dct < 1 - ct) {
		bt -= dbt;
		ct += dct;
		reti[0] = ai; reti[1] = bi; reti[2] = ci;
		rett[0] = bt; rett[1] = ct;
		return;
	    }

	    if ( (1 - ct)*dbt < bt*dct ) {
		t0 = (1 - ct)*dbt/dct;
		bx -= t0*ebx;
		by -= t0*eby;
		bt -= t0;
		ci += 1; 
		ct = 0;
		cx = polygon[(2*ci+0)%(2*sz)];
		cy = polygon[(2*ci+1)%(2*sz)];
		ecx = polygon[(2*ci+2)%(2*sz)] - polygon[(2*ci+0)%(2*sz)];
		ecy = polygon[(2*ci+3)%(2*sz)] - polygon[(2*ci+1)%(2*sz)];
	    } else {
		t0 = bt*dct/dbt;
		cx += t0*ecx;
		cy += t0*ecy;
		ct += t0;
		bi -= 1; 
		bt = 1;
		bx = polygon[(2*bi+2)%(2*sz)];
		by = polygon[(2*bi+3)%(2*sz)];
		ebx = polygon[(2*bi+2)%(2*sz)] - polygon[(2*bi+0)%(2*sz)];
		eby = polygon[(2*bi+3)%(2*sz)] - polygon[(2*bi+1)%(2*sz)];
	    }
	}
    }

    void maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status){
	//find a triangle inscribed in a convex polygon that is of maximal area among all inscribed triangles
	//input:
	//  polygon: a list of coordinates (x_0, y_0, x_1, y_1, ..., x_{n-1}, y_{n-1}) of the vertices of a convex polygon
	//           the vertices must be in couter-clockwise order and all must be extreme points
	//  ret[0] : is the index of the vertex to be used to obtain the initial anchored triangle
	//output:
	//  ret   : an array of 3 indices (i_a, i_b, i_c) between 0 and n-1
	//          the output triangle is (x_{i_a},y_{i_a})(x_{i_b},y_{i_b})(x_{i_c},y_{i_c})
	//  status : a flag indicating the successful termination of the function or an error encountered

	unsigned int ai,bi,ci;
	mpq_class bt, ct;
	unsigned int sz = polygon.size()/2;
	bool max_at_vert=false;
	unsigned int iter = 0;
	unsigned int maxiter = sz*10;
	unsigned int ai_start;
	mpq_class t_init[2];

	mpq_class ax,ay,bx,by,cx,cy;
	mpq_class eax,eay,ebx,eby,ecx,ecy;
	mpq_class tb1,tb2,tb3,tb4;
	mpq_class tc1,tc2,tc3,tc4;
	mpq_class qq, area;
	mpq_class amax = -1;

	ai = ret[0]%sz;
	eax = polygon[(2*ai+2)%(2*sz)] - polygon[(2*ai+0)%(2*sz)];
	eay = polygon[(2*ai+3)%(2*sz)] - polygon[(2*ai+1)%(2*sz)];
	anchored_triangle(polygon, -eay, eax, ret, t_init, status);
	ai = ret[0]; bi = ret[1]; ci = ret[2];
	bt = t_init[0]; ct = t_init[1];

	ai_start = ai;
	
	ax = polygon[(2*ai+0)%(2*sz)];
	ay = polygon[(2*ai+1)%(2*sz)];
	bx = polygon[(2*bi+0)%(2*sz)];
	by = polygon[(2*bi+1)%(2*sz)];
	cx = polygon[(2*ci+0)%(2*sz)];
	cy = polygon[(2*ci+1)%(2*sz)];

	eax = polygon[(2*ai+2)%(2*sz)] - polygon[(2*ai+0)%(2*sz)];
	eay = polygon[(2*ai+3)%(2*sz)] - polygon[(2*ai+1)%(2*sz)];
	ebx = polygon[(2*bi+2)%(2*sz)] - polygon[(2*bi+0)%(2*sz)];
	eby = polygon[(2*bi+3)%(2*sz)] - polygon[(2*bi+1)%(2*sz)];
	ecx = polygon[(2*ci+2)%(2*sz)] - polygon[(2*ci+0)%(2*sz)];
	ecy = polygon[(2*ci+3)%(2*sz)] - polygon[(2*ci+1)%(2*sz)];

	bx += bt*ebx;
	by += bt*eby;
	cx += ct*ecx;
	cy += ct*ecy;

	while (ai <= ai_start + sz) {

	    if (iter++ > maxiter) {*status=status_maxiter_exceeded; return; }

	    area = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
	    if ( (area > amax) || ( (area == amax) && !max_at_vert ) ) {
		amax = area;
		if (ct == 0 && bt == 0) {
		    max_at_vert = true;
		    ret[0] = (ai+sz)%sz;
		    ret[1] = (bi+sz)%sz;
		    ret[2] = (ci+sz)%sz;
		} else max_at_vert = false;
	    }
	    
	    if ((cy-by)*eax - (cx-bx)*eay == 0) { //support line coincides with forward edge at A, advance A
		ai += 1;
		ax = polygon[(2*ai+0)%(2*sz)];
		ay = polygon[(2*ai+1)%(2*sz)];
		eax = polygon[(2*ai+2)%(2*sz)] - polygon[(2*ai+0)%(2*sz)];
		eay = polygon[(2*ai+3)%(2*sz)] - polygon[(2*ai+1)%(2*sz)];
		continue;
	    }
	    if (bt > 1) {*status=status_runtime_error; return; }
	    if (ct > 1)  {*status=status_runtime_error; return; }
	    if (bt == 1) {
		bi += 1;
		bt = 0;
		bx = polygon[(2*bi+0)%(2*sz)];
		by = polygon[(2*bi+1)%(2*sz)];
		ebx = polygon[(2*bi+2)%(2*sz)] - polygon[(2*bi+0)%(2*sz)];
		eby = polygon[(2*bi+3)%(2*sz)] - polygon[(2*bi+1)%(2*sz)];
		continue;
	    }
	    if (ct == 1) {
		ci += 1;
		ct = 0;
		cx = polygon[(2*ci+0)%(2*sz)];
		cy = polygon[(2*ci+1)%(2*sz)];
		ecx = polygon[(2*ci+2)%(2*sz)] - polygon[(2*ci+0)%(2*sz)];
		ecy = polygon[(2*ci+3)%(2*sz)] - polygon[(2*ci+1)%(2*sz)];
		continue;
	    }

	    qq  = (cy-by)*((by-cy)*ebx*ecx - (ax-cx)*eby*ecx + (ax-bx)*ebx*ecy);
	    qq += (cx-bx)*((bx-cx)*eby*ecy - (ay-cy)*ebx*ecy + (ay-by)*eby*ecx);
	    
	    if (qq > 0) {
		//C stays fixed, B moves along eb
		//possible stopping conditions: 
		//   (1) B hits end of edge,
		//   (3) BC becomes parallel to ea,
		//   (4) Q becomes 0
		tb1 = 1 - bt;

		if (eay*ebx - eax*eby == 0) tb3 = -1;
		else tb3 = (eax*(cy-by) - eay*(cx-bx))/(eax*eby - eay*ebx);

		if (ebx*ecy - ecx*eby == 0) tb4 = -1;
		else {
		    tb4  = ebx*ecx*(by-cy)*(by-cy);
		    tb4 += eby*ecy*(bx-cx)*(bx-cx);
		    tb4 += ebx*ecy*( (ax-bx)*by + (cx-bx)*ay + (2*bx-ax-cx)*cy);
		    tb4 += eby*ecx*( (bx-cx)*ay - (cx-ax)*cy + (2*cx-ax-bx)*by);
		    tb4 /= ebx*ecy - eby*ecx;
		    tb4 /= ebx*(ay + by - 2*cy) - eby*(ax + bx - 2*cx);
		}

		if (tb4 < tb1 && tb4 > 0) tb1 = tb4;
		if (tb3 < tb1 && tb3 > 0) tb1 = tb3;

		if (tb1 < 0) {*status=status_runtime_error; return; }

		bx += tb1*ebx;
		by += tb1*eby;
		bt += tb1;

		continue;
	    }

	    if (qq < 0) {
		//B stays fixed, C moves along ec
		//possible stopping conditions: 
		//   (2) C hits end of edge,
		//   (3) BC becomes parallel to ea,
		//   (4) Q becomes 0
		tc2 = 1 - ct;

		if (eay*ecx - eax*ecy == 0) tc3 = -1;
		else tc3 = (eax*(by-cy) - eay*(bx-cx))/(eax*ecy - eay*ecx);

		if (ebx*ecy - ecx*eby == 0) tc4 = -1;
		else {
		    tc4  = ebx*ecx*(by-cy)*(by-cy);
		    tc4 += eby*ecy*(bx-cx)*(bx-cx);
		    tc4 += ebx*ecy*( (ax-bx)*by + (cx-bx)*ay + (2*bx-ax-cx)*cy);
		    tc4 += eby*ecx*( (bx-cx)*ay - (cx-ax)*cy + (2*cx-ax-bx)*by);
		    tc4 /= -( ebx*ecy - eby*ecx );
		    tc4 /= ecx*(ay + cy - 2*by) - ecy*(ax + cx - 2*bx);
		}

		if (tc4 < tc2 && tc4 > 0) tc2 = tc4;
		if (tc3 < tc2 && tc3 > 0) tc2 = tc3;
		
		if (tc2 < 0) {*status=status_runtime_error; return; }

		cx += tc2*ecx;
		cy += tc2*ecy;
		ct += tc2;

		continue;
	    }

	    if (qq == 0) {
		//B and C both move
		//possible stopping conditions: 
		//   (1) B hits end of edge,
		//   (2) C hits end of edge,
		//   (3) BC becomes parallel to ea,

		mpq_class beta =  ebx*(2*cy-ay-by) - eby*(2*cx-bx-ax);
		mpq_class gamma = ecx*(2*by-ay-cy) - ecy*(2*bx-cx-ax);
		mpq_class alpha = 2*(ebx*ecy - ecx*eby);

		tb1 = 1 - bt;
		if (gamma - alpha*tb1 == 0) {tb1 = tc1 = -1;}
		else tc1 = beta*tb1/(gamma - alpha*tb1);

		tc2 = 1 - ct;
		if (beta + alpha*tc2 == 0) {tb2 = tc2 = -1;}
		else tb2 = gamma*tc2/(beta + alpha*tc2);

		if (eay*ebx - eax*eby == 0) {tb3 = tc3 = -1;}
		else if (eay*ecx - eax*ecy == 0) {tb3 = tc3 = -1;}
		else {
		    tb3  = ebx*ecx*eay*(by-cy)+eby*ecy*eax*(bx-cx);
		    tb3 += ebx*ecy*(eax*(by-ay) - eay*(2*bx-ax-cx));
		    tb3 += eby*ecx*(eay*(bx-ax) - eax*(2*by-ay-cy));
		    tb3 /= alpha*(eay*ebx - eax*eby);

		    tc3  = ebx*ecx*eay*(by-cy)+eby*ecy*eax*(bx-cx);
		    tc3 += eby*ecx*(-eax*(cy-ay) + eay*(2*cx-ax-bx));
		    tc3 += ebx*ecy*(-eay*(cx-ax) + eax*(2*cy-ay-by));
		    tc3 /= alpha*(eay*ecx - eax*ecy);
		}

		if ( (tc1 < 0) || (tb2 <= tb1 && tc2 <= tc1) ) { tc1 = tc2; tb1 = tb2;}
		if (tb3 <= tb1 && tc3 <= tc1 && tb3 >= 0 && tc3 >= 0) {tb1 = tb3; tc1 = tc3;}
		if (tb1 < 0 || tc1 < 0) {*status=status_runtime_error; return; }

		bx += tb1*ebx;
		by += tb1*eby;
		bt += tb1;
		cx += tc1*ecx;
		cy += tc1*ecy;
		ct += tc1;
	    }
	}

	if (!max_at_vert) *status = status_max_not_at_vertex;
	return; 
    }

    void brute_force_maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status){
	mpq_class amax, area;
	mpq_class ax,ay,bx,by,cx,cy;
	unsigned int sz = polygon.size()/2;
	unsigned int ai,bi,ci;

	amax = -1;
	*status = status_ok;
	for (ai = 0; ai+2 < sz; ai++) { 
	    ax = polygon[(2*ai+0)%(2*sz)];
	    ay = polygon[(2*ai+1)%(2*sz)];
	    for (bi = ai+1; bi+1 < sz; bi++) { 
		bx = polygon[(2*bi+0)%(2*sz)];
		by = polygon[(2*bi+1)%(2*sz)];
		for (ci = bi+1; ci < sz; ci++) { 
		    cx = polygon[(2*ci+0)%(2*sz)];
		    cy = polygon[(2*ci+1)%(2*sz)];
		    area = (bx-ax)*(cy-ay) - (cx-ax)*(by-ay);
		    if (area > amax) {
			amax = area;
			ret[0] = ai;
			ret[1] = bi;
			ret[2] = ci;
		    }
		}
	    }
	}
	return;
    }
}
