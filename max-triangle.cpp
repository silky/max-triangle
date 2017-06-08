#include <iomanip>
#include "max-triangle.hpp"

namespace InscribedTriangle {
    void print_status(state_flag status) {
	switch(status) {
	    case status_ok: std::cout << "status_ok" << std::endl; break;
	    case status_not_convex: std::cout << "status_not_convex" << std::endl; break;
	    case status_no_interior: std::cout << "status_no_interior" << std::endl; break;
	    case status_runtime_error: std::cout << "status_runtime_error" << std::endl; break;
	    case status_max_not_at_vertex: std::cout << "status_max_not_at_vertex" << std::endl; break;
	    case status_maxiter_exceeded: std::cout << "status_maxiter_exceeded" << std::endl; break;
	}
    }

    bool is_convex(std::vector<mpq_class> polygon) {
	//check if polygon is a strictly convex polygon with vertices listed in CCW order

	unsigned int i;
	unsigned int sz = polygon.size()/2;
	mpq_class x;

	for (i = 0; i<polygon.size(); i++) {
	    x  = (polygon[(2*i+2)%(2*sz)] - polygon[(2*i+0)%(2*sz)])*(polygon[(2*i+5)%(2*sz)] - polygon[(2*i+1)%(2*sz)]);
	    x -= (polygon[(2*i+4)%(2*sz)] - polygon[(2*i+0)%(2*sz)])*(polygon[(2*i+3)%(2*sz)] - polygon[(2*i+1)%(2*sz)]);
	    if (x <= 0) return false;
	}

	return true;

    }


    void anchored_triangle(std::vector<mpq_class> polygon, mpq_class nx, mpq_class ny, unsigned int *reti, mpq_class *rett, state_flag *status){

	mpq_class minx,maxx,x;
	unsigned int ai,bi,ci;
	mpq_class ax,ay,bx,by,cx,cy;
	mpq_class ebx,eby,ecx,ecy;
	mpq_class dbt, dct, t0, tq;
	mpq_class ee, neb, nec;
	mpq_class bt,ct;
	unsigned int i;
	unsigned int sz = polygon.size()/2;

	if (polygon.size() < 3) { *status = status_no_interior; return; }
	if (!is_convex(polygon)) { *status = status_not_convex; return; }
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

	    //dbt  = nx*eby*ecx*(ax-bx) + nx*ecy*ebx*(2*bx-ax-cx);
	    //dbt += ny*ebx*ecy*(by-ay) + ny*ecx*eby*(ay+cy-2*by);
	    //dbt /= 2*(ebx*ecy - ecx*eby)*(ebx*nx + eby*ny);

	    //dct  = nx*ecy*ebx*(ax-cx) + nx*eby*ecx*(2*cx-ax-bx);
	    //dct += ny*ecx*eby*(cy-ay) + ny*ebx*ecy*(ay+by-2*cy);
	    //dct /= 2*(ebx*ecy - ecx*eby)*(ecx*nx + ecy*ny);
	    
	    tq  = nx*( (by-cy)*ebx*ecx + (cx-ax)*eby*ecx + (ax-bx)*ebx*ecy);
	    tq -= ny*( (bx-cx)*eby*ecy + (cy-ay)*ebx*ecy + (ay-by)*eby*ecx);
	    //tq -= ny*(ay*eby*ecx - by*eby*ecx - ay*ebx*ecy + cy*ebx*ecy + bx*eby*ecy - cx*eby*ecy);
	    tq /= 2*ee;

	    dbt = -tq/neb;
	    dct = tq/nec;

	    //std::cout << bi << ", " << ci << ", " << bt << ", " << ct << ", " << dbt << ", " << dct << ", " << bx << ", " << by << ", " << cx << ", " << cy << std::endl;

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

	//mpq_class nx, ny;
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

	/*if (polygon.size() < 3) { *status = status_no_interior; return;  }
	if (!is_convex(polygon)) { *status = status_not_convex; return;  }
	*status = status_ok;

	bi = ret[0];
	ci = bi+1;
	ai = bi+2;
	bt = 0;
	ct = 0;
	area = (polygon[(2*ci+0)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ai+1)%(2*sz)]-polygon[(2*bi+1)%(2*sz)]) - (polygon[(2*ai+0)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ci+1)%(2*sz)]-polygon[(2*bi+1)%(2*sz)]);
	amax = area;
	while (true){
	    if (area <= (polygon[(2*ci+2)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ai+1)%(2*sz)]-polygon[(2*bi+1)%(2*sz)]) - (polygon[(2*ai+0)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ci+3)%(2*sz)]-polygon[(2*bi+1)%(2*sz)])) ci++;
	    else if (area <= (polygon[(2*ci+0)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ai+3)%(2*sz)]-polygon[(2*bi+1)%(2*sz)]) - (polygon[(2*ai+2)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ci+1)%(2*sz)]-polygon[(2*bi+1)%(2*sz)])) ai++;
	    else break;
	    area = (polygon[(2*ci+0)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ai+1)%(2*sz)]-polygon[(2*bi+1)%(2*sz)]) - (polygon[(2*ai+0)%(2*sz)]-polygon[(2*bi+0)%(2*sz)])*(polygon[(2*ci+1)%(2*sz)]-polygon[(2*bi+1)%(2*sz)]);
	}*/

	ai = ret[0];
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

	//nx = cy - by;
	//ny = bx - cx;

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

#ifdef DEBUG
	    std::cout << ai << "\t" << bi << "\t" << ci << "\t";
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << bt.get_d() << "\t";
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << ct.get_d() << "\t";
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << area.get_d() << std::endl;
#endif


	    //ee = ebx*ecy - ecx*eby;
	    //fb = ebx*(ay + by - 2*cy) - eby*(ax + bx - 2*cx);
	    //fc = ecx*(ay + cy - 2*by) - ecy*(ax + cx - 2*bx);
	    //tq  = nx*(by*ebx*ecx - cy*ebx*ecx - ax*eby*ecx + cx*eby*ecx + ax*ebx*ecy - bx*ebx*ecy);
	    //tq -= ny*(ay*eby*ecx - by*eby*ecx - ay*ebx*ecy + cy*ebx*ecy + bx*eby*ecy - cx*eby*ecy);
	    qq  = (cy-by)*((by-cy)*ebx*ecx - (ax-cx)*eby*ecx + (ax-bx)*ebx*ecy);
	    qq += (cx-bx)*((bx-cx)*eby*ecy - (ay-cy)*ebx*ecy + (ay-by)*eby*ecx);
#ifdef DEBUG
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << qq.get_d() << "\t";
#endif
	    //qq > 0 ---> move right
	    //qq < 0 ---> move left
	    //qq = 0 ---> swing
	    
	    if (qq > 0) {
		//C stays fixed, B moves along eb
		//possible stopping conditions: 
		//   (1) B hits end of edge,
		//   (3) BC becomes parallel to ea,
		//   (4) Q becomes 0
		tb1 = 1 - bt;

		if (eay*ebx - eax*eby == 0) tb3 = -1;
		//else t2 = (by*eax - cy*eax - bx*eay + cx*eay)/(eay*ebx - eax*eby);
		else tb3 = (eax*(cy-by) - eay*(cx-bx))/(eax*eby - eay*ebx);

		if (ebx*ecy - ecx*eby == 0) tb4 = -1;
		else {
		    tb4  = ebx*ecx*(by-cy)*(by-cy);
		    tb4 += eby*ecy*(bx-cx)*(bx-cx);
		    //t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
		    //t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
		    tb4 += ebx*ecy*( (ax-bx)*by + (cx-bx)*ay + (2*bx-ax-cx)*cy);
		    tb4 += eby*ecx*( (bx-cx)*ay - (cx-ax)*cy + (2*cx-ax-bx)*by);
		    tb4 /= ebx*ecy - eby*ecx;
		    tb4 /= ebx*(ay + by - 2*cy) - eby*(ax + bx - 2*cx);
		}

#ifdef DEBUG
		std::cout << "qq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qq.get_d() << "\t";
		std::cout << "t1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1.get_d() << "\t";
		std::cout << "t2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2.get_d() << "\t";
		std::cout << "t3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3.get_d() << std::endl;
#endif

		if (tb4 < tb1 && tb4 > 0) tb1 = tb4;
		if (tb3 < tb1 && tb3 > 0) tb1 = tb3;

		if (tb1 < 0) {*status=status_runtime_error; return; }

		bx += tb1*ebx;
		by += tb1*eby;
		bt += tb1;
		//nx = cy - by;
		//ny = bx - cx;

		continue;
	    }

	    if (qq < 0) {
		//B stays fixed, C moves along ec
		//possible stopping conditions: 
		//   (2) C hits end of edge,
		//   (3) BC becomes parallel to ea,
		//   (4) Q becomes 0
		tc2 = 1 - ct;

		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << ecy.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << eax.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << eay.get_d() << std::endl;
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << cx.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << cy.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << bx.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << by.get_d() << std::endl;
		if (eay*ecx - eax*ecy == 0) tc3 = -1;
		//else t2 = -(by*eax - cy*eax - bx*eay + cx*eay)/(eay*ecx - eax*ecy);
		else tc3 = (eax*(by-cy) - eay*(bx-cx))/(eax*ecy - eay*ecx);

		if (ebx*ecy - ecx*eby == 0) tc4 = -1;
		else {
		    tc4  = ebx*ecx*(by-cy)*(by-cy);
		    tc4 += eby*ecy*(bx-cx)*(bx-cx);
		    //t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
		    //t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
		    tc4 += ebx*ecy*( (ax-bx)*by + (cx-bx)*ay + (2*bx-ax-cx)*cy);
		    tc4 += eby*ecx*( (bx-cx)*ay - (cx-ax)*cy + (2*cx-ax-bx)*by);
		    tc4 /= -( ebx*ecy - eby*ecx );
		    tc4 /= ecx*(ay + cy - 2*by) - ecy*(ax + cx - 2*bx);
		}

#ifdef DEBUG
		std::cout << "qq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qq.get_d() << "\t";
		std::cout << "t1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1.get_d() << "\t";
		std::cout << "t2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2.get_d() << "\t";
		std::cout << "t3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3.get_d() << std::endl;
#endif


		if (tc4 < tc2 && tc4 > 0) tc2 = tc4;
		if (tc3 < tc2 && tc3 > 0) tc2 = tc3;
		
		if (tc2 < 0) {*status=status_runtime_error; return; }

		cx += tc2*ecx;
		cy += tc2*ecy;
		ct += tc2;
		//nx = cy - by;
		//ny = bx - cx;

		continue;
	    }

	    if (qq == 0) {
		//B and C both move
		//possible stopping conditions: 
		//   (1) B hits end of edge,
		//   (2) C hits end of edge,
		//   (3) BC becomes parallel to ea,
#ifdef DEBUG
		mpq_class qu = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby);
		std::cout << "alpha1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qu.get_d() << "\t";
		qu = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy);
		std::cout << "alpha2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qu.get_d() << "\t";
		qu = (eay*ebx - eax*eby);
		std::cout << "eaxab: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qu.get_d() << "\t";
		std::cout << "ee: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << ee.get_d() << std::endl;
#endif


		//mpq_class beta = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby);
		//mpq_class gamma = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy);
		mpq_class beta =  ebx*(2*cy-ay-by) - eby*(2*cx-bx-ax);
		mpq_class gamma = ecx*(2*by-ay-cy) - ecy*(2*bx-cx-ax);
		mpq_class alpha = 2*(ebx*ecy - ecx*eby);

		//time in terms of b's motion until it hits vertex
		tb1 = 1 - bt;
		//if (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*tb1 == 0) {tb1 = tc1 = -1;}
		if (gamma - alpha*tb1 == 0) {tb1 = tc1 = -1;}
		else tc1 = beta*tb1/(gamma - alpha*tb1);

		//time in terms of c's motion until it hits vertex
		tc2 = 1 - ct;
		//time in terms of b's motion until c hits vertex
		//if (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*tc2 == 0) {tb2 = tc2 = -1;}
		//else tb2 = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy)*tc2/(ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*tc2);
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

#ifdef DEBUG
		std::cout << "qq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qq.get_d() << "\t";
		std::cout << "tb1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tb1.get_d() << "\t";
		std::cout << "tb2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tb2.get_d() << "\t";
		std::cout << "tb3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tb3.get_d() << std::endl;
		std::cout << "qq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << qq.get_d() << "\t";
		std::cout << "tc1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tc1.get_d() << "\t";
		std::cout << "tc2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tc2.get_d() << "\t";
		std::cout << "tc3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tc3.get_d() << std::endl;
#endif

		if ( (tc1 < 0) || (tb2 <= tb1 && tc2 <= tc1) ) { tc1 = tc2; tb1 = tb2;}
		if (tb3 <= tb1 && tc3 <= tc1 && tb3 >= 0 && tc3 >= 0) {tb1 = tb3; tc1 = tc3;}
		if (tb1 < 0 || tc1 < 0) {*status=status_runtime_error; return; }

		bx += tb1*ebx;
		by += tb1*eby;
		bt += tb1;
		cx += tc1*ecx;
		cy += tc1*ecy;
		ct += tc1;
		//nx = cy - by;
		//ny = bx - cx;
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
