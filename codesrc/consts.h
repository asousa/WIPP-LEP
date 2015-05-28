
// Constants used throughout simulation

#define		PI		3.14159265358979311599796346854418516159
#define		D2R		PI/180.0
#define		R2D		180.0/PI

#define		NUM_STEPS	15000
#define		T_STEP		0.0004


#define		Q_EL		1.602E-19
#define		M_EL		9.1E-31
#define		E_EL		5.105396765648739E5
#define		MU0		PI*4E-7		
#define		EPS0		8.854E-12
#define		C		2.997956376932163e+08
#define		Z0		377.0
#define		R_E		6370000.0
#define		H_MAGNETO	1E6
#define		H_IONO		1E5

#define		A		5E3	
#define		B		1E5
#define		H_E		5000.0
#define		I0		-100000.0

#define		P_DIST		0.0
#define		Q_DIST		2.0
#define		AN_CM_DIST	2E5
#define		V0_DIST		1.0
#define		M_RES		0

#define		RES_DT		0.02
#define		RES_FINT	5.0
#define		C		2.997956376932163e+08
#define		R_E		6370000.0
#define		H_IONO		1E5
#define		M_EL		9.1E-31
#define		E_EL		5.105396765648739E5

#define		DE_EXP		0.003
#define		E_EXP_BOT	1.477
#define		E_EXP_TOP	7.477
#define		NUM_E		2200
#define		SQUARE		1

#define		EALimS		-40.0
#define		EALimN		40.0
#define		EAIncr		1.0
#define		dL0		6E-4
#define		DF		50.0
#define		DT		0.01
#define		NUMLATS		((EALimN - EALimS)/EAIncr + 1)

#define		EA_SPLIT	1
#define		MULT		2.0


long freqs[] = {	200,
      210,	
			220,		
			230,
			240,
			250,
			260,
			270,
			280,
			290,
			
			300,
			310, 
		       	320,
			330,
			340,
			350,
			360,
			370,
			380,
			390,
			
			400,
			410,
			420,
			430,
			440,
			460,
			480,
			500,
			520,
			540,
			
			560,
			580,
			600,
			620,
			640,
			660,
			690,
			720,
			750,
			780, 
			
			810,
			840,
			880,
			920,
			960, 
			1000,  
			1040,
			1080,
			1120,
			1160,
			
			1200,
			1250,
			1300,
			1360,
			1420,
			1480,
			1540,
			1600,
			1670,
			1740,
				
			1800,
			1880,
			1960,
			2050,
			2140,
			2230,
			2320,
			2410,
			2500,
			2600,
			
			2700,
			2800,
			2900, 
			3000,
			3100,
			3200,
			3300,
			3450,
			3600,
			3750,
			
			3900,
			4050,
			4200,
			4400,
			4600, 
			4800,
			5000,
			5200,
			5400,
			5600,
			
			5800,
			6000,
			6200,
			6400,
			6600,
			6900,
			7200,
			7500,
			7800,
			8100,
			
			8400,
			8800,
			9200,
			9600,
			10000,
			10500,
			11000,
			11500,
			12000,
			12500,
			
			13000,
			13500,
			14000,
			14500,
			15000,
			16000,
			17000,
			18000,
			19000,
			20000,
			
			22000,
			24000,
			26000,
			28000,
			30000,
			35000,
			40000,
			45000,
			50000,
			60000     };

// Freqs: 32, log-spaced between 200 and 60000
// long freqs[] = { 200,
// 			240,
// 			289,
// 			347,
// 			418,
// 			502,
// 			603,
// 			725,
// 			872,
// 			1048,
// 			1259,
// 			1514,
// 			1819,
// 			2187,
// 			2629,
// 			3160,
// 			3798,
// 			4565,
// 			5487,
// 			6596,
// 			7928,
// 			9530,
// 			11455,
// 			13769,
// 			16550,
// 			19893,
// 			23912,
// 			28742,
// 			34549,
// 			41528,
// 			49916,
// 			60000 };





// Definitions: (damping.c)
//-------------

// #define		PI		3.14159265358979311599796346854418516159
// #define		D2R		PI/180.0
// #define		R2D		180.0/PI

// #define		NUM_STEPS	15000
// #define		T_STEP		0.0004

// #define		Q_EL		1.602E-19
// #define		M_EL		9.1E-31
// #define		MU0		PI*4E-7		
// #define		EPS0		8.854E-12
// #define		C		2.997956376932163e+08
// #define		Z0		377.0
// #define		R_E		6370000.0
// #define		H_MAGNETO	1E6
// #define		H_IONO		1E5

// #define		A		5E3	
// #define		B		1E5
// #define		H_E		5000.0
// #define		I0		-10530.0

// #define		P_DIST		0.0
// #define		Q_DIST		2.0
// #define		AN_CM_DIST	2E5
// #define		V0_DIST		1.0
// #define		M_RES		0





// DEFINITIONS (calcFlux_austin.d)
// -----------


// #define		RES_DT		0.02
// #define		RES_FINT	5.0
// #define		C		2.997956376932163e+08
// #define		R_E		6370000.0
// #define		H_IONO		1E5
// #define		PI		3.141592653589793115997963468544185
// #define		D2R		PI/180.0
// #define		M_EL		9.1E-31
// #define		E_EL		5.105396765648739E5

// #define		DE_EXP		0.003
// #define		E_EXP_BOT	1.477
// #define		E_EXP_TOP	7.477
// #define		NUM_E		2200
// #define		SQUARE		1



// Definitions: (readone.c)
//-------------

// #define		PI		3.141592653589793115997963468544185161590576
// #define		D2R		PI/180.0
// #define		R2D		180.0/PI

// #define		NUM_STEPS	15000
// #define		T_STEP		0.0004

// #define		Q_EL		1.602E-19
// #define		M_EL		9.1E-31
// #define		E_EL		5.105396765648739E5
// #define		MU0		PI*4E-7		
// #define		EPS0		8.854E-12
// #define		C		2.997956376932163e+08
// #define		Z0		377.0
// #define		R_E		6370000.0
// #define		H_MAGNETO	1E6
// #define		H_IONO		1E5

// #define		A		5E3	
// #define		B		1E5
// #define		H_E		5000.0
// #define		I0		-100000.0
// // peak current defined here, in Amps!

// #define		P_DIST		0.0
// #define		Q_DIST		2.0
// #define		AN_CM_DIST	2E5
// #define		V0_DIST		1.0
// #define		M_RES		0

// #define		EALimS		-40.0
// #define		EALimN		40.0
// #define		EAIncr		1.0
// #define		dL0		6E-4
// #define		DF		50.0
// #define		DT		0.01
// #define		NUMLATS		((EALimN - EALimS)/EAIncr + 1)

// // #define		DV_TOT		1e5
// #define		DE_EXP		0.003
// #define		E_EXP_BOT	1.477
// #define		E_EXP_TOP	7.477	
// #define		RES_DT		0.02
// #define		RES_FINT	5.0
// #define		EA_SPLIT	1
// #define		MULT		2.0
// #define		NUM_E		2200
