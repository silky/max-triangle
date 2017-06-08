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

    void maximum_triangle(std::vector<mpq_class> polygon, unsigned int *ret, state_flag *status){

	mpq_class nx, ny;
	unsigned int ai,bi,ci;
	mpq_class bt, ct;
	unsigned int sz = polygon.size()/2;
	bool max_at_vert=false;
	unsigned int iter = 0;
	unsigned int maxiter = sz*10;
	unsigned int ai_start;

	mpq_class ax,ay,bx,by,cx,cy;
	mpq_class pax,pay,pbx,pby,pcx,pcy;
	mpq_class eax,eay,ebx,eby,ecx,ecy;
	mpq_class t1,t2,t3;
	mpq_class ee,fb,fc,tq;
	mpq_class area;
	mpq_class amax = -1;

	if (polygon.size() < 3) { *status = status_no_interior; return;  }
	if (!is_convex(polygon)) { *status = status_not_convex; return;  }
	*status = status_ok;

	bi = 0;
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
	}

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

	nx = cy - by;
	ny = bx - cx;

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
	    
	    if (nx*eax + ny*eay == 0) { //support line coincides with forward edge at A, advance A
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


	    ee = ebx*ecy - ecx*eby;
	    fb = ebx*(ay + by - 2*cy) - eby*(ax + bx - 2*cx);
	    fc = ecx*(ay + cy - 2*by) - ecy*(ax + cx - 2*bx);
	    tq  = nx*(by*ebx*ecx - cy*ebx*ecx - ax*eby*ecx + cx*eby*ecx + ax*ebx*ecy - bx*ebx*ecy);
	    tq -= ny*(ay*eby*ecx - by*eby*ecx - ay*ebx*ecy + cy*ebx*ecy + bx*eby*ecy - cx*eby*ecy);
#ifdef DEBUG
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << ee.get_d() << "\t";
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << fb.get_d() << "\t";
	    std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << fc.get_d() << std::endl;
#endif
	    //tq > 0 ---> move right
	    //tq < 0 ---> move left
	    //tq = 0 ---> swing
	    
	    if (tq > 0) {
		//C stays fixed, B moves along eb
		//possible stopping conditions: 
		//   (1) B hits end of edge,
		//   (2) BC becomes parallel to ea,
		//   (3) tq becomes 0
		t1 = 1 - bt;

		if (eay*ebx - eax*eby == 0) t2 = -1;
		else t2 = (by*eax - cy*eax - bx*eay + cx*eay)/(eay*ebx - eax*eby);

		if (ee == 0) t3 = -1;
		else {
		    t3  = ebx*ecx*(by-cy)*(by-cy);
		    t3 += eby*ecy*(bx-cx)*(bx-cx);
		    t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
		    t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
		    t3 /= ee*fb;
		}

#ifdef DEBUG
		std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
		std::cout << "t1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1.get_d() << "\t";
		std::cout << "t2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2.get_d() << "\t";
		std::cout << "t3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3.get_d() << std::endl;
#endif

		if (t3 < t1 && t3 > 0) t1 = t3;
		if (t2 < t1 && t2 > 0) t1 = t2;

		if (t1 < 0) {*status=status_runtime_error; return; }

		bx += t1*ebx;
		by += t1*eby;
		bt += t1;
		nx = cy - by;
		ny = bx - cx;

		continue;
	    }

	    if (tq < 0) {
		//B stays fixed, C moves along ec
		//possible stopping conditions: 
		//   (1) C hits end of edge,
		//   (2) BC becomes parallel to ea,
		//   (3) tq becomes 0
		t1 = 1 - ct;

		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << ecy.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << eax.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << eay.get_d() << std::endl;
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << cx.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << cy.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << bx.get_d() << "\t";
		//std::cout << std::fixed << std::setfill(' ') << std::setw(10) << std::setprecision(3) << by.get_d() << std::endl;
		if (eay*ecx - eax*ecy == 0) t2 = -1;
		else t2 = -(by*eax - cy*eax - bx*eay + cx*eay)/(eay*ecx - eax*ecy);

		if (ee == 0) t3 = -1;
		else {
		    t3  = ebx*ecx*(by-cy)*(by-cy);
		    t3 += eby*ecy*(bx-cx)*(bx-cx);
		    t3 += ebx*ecy*(-ay*bx + ax*by - bx*by + ay*cx - ax*cy + 2*bx*cy - cx*cy);
		    t3 += eby*ecx*( ay*bx - ax*by - bx*by - ay*cx + 2*by*cx + ax*cy - cx*cy);
		    t3 /= -ee*fc;
		}

#ifdef DEBUG
		std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
		std::cout << "t1: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1.get_d() << "\t";
		std::cout << "t2: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2.get_d() << "\t";
		std::cout << "t3: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3.get_d() << std::endl;
#endif


		if (t3 < t1 && t3 > 0) t1 = t3;
		if (t2 < t1 && t2 > 0) t1 = t2;
		
		if (t1 < 0) {*status=status_runtime_error; return; }

		cx += t1*ecx;
		cy += t1*ecy;
		ct += t1;
		nx = cy - by;
		ny = bx - cx;

		continue;
	    }

	    if (tq == 0) {
		mpq_class t1b,t1c,t2b,t2c,t3b,t3c;
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

		//time in terms of b's motion until it hits vertex
		t1b = 1 - bt;
		if (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*t1b == 0) {t1b = t1c = -1;}
		else t1c = (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby)*t1b/(ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy + 2*ee*t1b);

		//time in terms of c's motion until it hits vertex
		t2c = 1 - ct;
		//time in terms of b's motion until c hits vertex
		if (ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*t2c == 0) {t2b = t2c = -1;}
		else t2b = (ay*ecx - 2*by*ecx + cy*ecx - ax*ecy + 2*bx*ecy - cx*ecy)*t2c/(ay*ebx + by*ebx - 2*cy*ebx - ax*eby - bx*eby + 2*cx*eby - 2*ee*t2c);

		if (eay*ebx - eax*eby == 0) {t3b = t3c = -1;}
		else if (eay*ecx - eax*ecy == 0) {t3b = t3c = -1;}
		else {
		    t3b  = ebx*ecx*eay*(by-cy);
		    t3b += ebx*ecy*(-ay*eax + by*eax + ax*eay - 2*bx*eay + cx*eay);
		    t3b += eby*ecx*( ay*eax - 2*by*eax + cy*eax - ax*eay + bx*eay);
		    t3b += eby*ecy*eax*(bx-cx);
		    t3b /= 2*ee*(eay*ebx - eax*eby);

		    t3c  = ebx*ecx*eay*(by-cy);
		    t3c += ebx*ecy*(-ay*eax - by*eax + 2*cy*eax + ax*eay - cx*eay);
		    t3c += eby*ecx*( ay*eax - cy*eax - ax*eay - bx*eay + 2*cx*eay);
		    t3c += eby*ecy*eax*(bx-cx);
		    t3c /= 2*ee*(eay*ecx - eax*ecy);
		}

#ifdef DEBUG
		std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
		std::cout << "t1b: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1b.get_d() << "\t";
		std::cout << "t2b: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2b.get_d() << "\t";
		std::cout << "t3b: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3b.get_d() << std::endl;
		std::cout << "tq: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << tq.get_d() << "\t";
		std::cout << "t1c: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t1c.get_d() << "\t";
		std::cout << "t2c: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t2c.get_d() << "\t";
		std::cout << "t3c: " << std::scientific << std::setfill(' ') << std::scientific << std::setw(10) << std::setprecision(3) << t3c.get_d() << std::endl;
#endif

		if ( (t1c < 0) || (t2b <= t1b && t2c <= t1c) ) { t1c = t2c; t1b = t2b;}
		if (t3b <= t1b && t3c <= t1c && t3b >= 0 && t3c >= 0) {t1b = t3b; t1c = t3c;}
		if (t1b < 0 || t1c < 0) {*status=status_runtime_error; return; }

		bx += t1b*ebx;
		by += t1b*eby;
		bt += t1b;
		cx += t1c*ecx;
		cy += t1c*ecy;
		ct += t1c;
		nx = cy - by;
		ny = bx - cx;
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
