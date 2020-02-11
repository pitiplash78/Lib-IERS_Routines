using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IERS_Routines
{
    public static class CommonFunctions
    {
        public static void fundarg(double T, ref double L, ref double LP, ref double F, ref double D, ref double OM) // SUBROUTINE FUNDARG ( T, L, LP, F, D, OM )
        {
            /*   F U N D A R G
            
             * *  This routine is part of the International Earth Rotation and
             *  Reference Systems Service (IERS) Conventions software collection.

             *  This subroutine computes the lunisolar fundamental arguments.
             *  The model used is from Simon et al. (1994) as recommended by the IERS
             *  Conventions (2010).  Refer to IERS Conventions (2010) Chapter 5
             *  Sections 5.7.1 - 5.7.2 (pp. 57 - 59).

             *  In general, Class 1, 2, and 3 models represent physical effects that
             *  act on geodetic parameters while canonical models provide lower-level
             *  representations or basic computations that are used by Class 1, 2, or
             *  3 models.

             *  Status: Canonical model

             *     Class 1 models are those recommended to be used a priori in the
             *     reduction of raw space geodetic data in order to determine
             *     geodetic parameter estimates.
             *     Class 2 models are those that eliminate an observational
             *     singularity and are purely conventional in nature.
             *     Class 3 models are those that are not required as either Class
             *     1 or 2.
             *     Canonical models are accepted as is and cannot be classified as a
             *     Class 1, 2, or 3 model.

             *  Given:
             *     T           d      TT, Julian centuries since J2000 (Note 1)

             *  Returned:
             *     L           d      Mean anomaly of the Moon (Note 2)
             *     LP          d      Mean anomaly of the Sun (Note 2)
             *     F           d      L - OM (Notes 2 and 3
             *     D           d      Mean elongation of the Moon from the Sun (Note 2) 
             *     OM          d      Mean longitude of the ascending node of the Moon (Note 2)

             *  Notes:
             *  1) Though T is strictly TDB, it is usually more convenient to use TT, which makes no significant difference. Julian centuries since 
             *     J2000 is (JD - 2451545.0)/36525.

             *  2) The expression used is as adopted in IERS Conventions (2010) and is from Simon et al. (1994).  Arguments are in radians.

             *  3) L in this instance is the Mean Longitude of the Moon. OM is the Mean longitude of the ascending node of the Moon.

             *  Test case:
             *     given input: T = 0.07995893223819302 Julian centuries since J2000
             *                  (MJD = 54465)
             *     expected output:  L = 2.291187512612069099 radians
             *                       LP = 6.212931111003726414 radians
             *                       F = 3.658025792050572989 radians
             *                       D = 4.554139562402433228 radians
             *                       OM = -0.5167379217231804489 radians

             *  References:
             *     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M., Francou, G., Laskar, J., 1994, Astron.Astrophys. 282, 663-683
             *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010), IERS Technical Note No. 36, BKG (2010)

             * Revisions: 
             * 2008 January 18 B.E.Stetzler  Initial changes to header and used 2PI instead of PI as parameter
             * 2008 January 25 B.E. Stetzler Additional changes to header and defined fundamental arguments
             * 2008 January 28 B.E. Stetzler Additional changes to header 
             * 2008 March   12 B.E. Stetzler Applied changes to wording of notes.
             * 2008 April   03 B.E. Stetzler Provided example test case.
             * 2009 February 11 B.E. Stetzler Corrected term in OM from 6962890.2665 to 6962890.5431 and updated test case
             * 2009 May     07 B.E. Stetzler Code formatting changes based on client recommendations
             * 2009 May     07 B.E. Stetzler Updated test case due to above changes 
             * 2010 February 25 B.E. Stetzler Recalculation of fundamental arguments */



            /*  Compute the fundamental argument L. */
            /*    L = MOD (       485868.249036D0 +
                 .                  T*( 1717915923.2178D0 +
                 .                  T*(         31.8792D0 +
                 .                  T*(          0.051635D0 +
                 .                  T*(        - 0.00024470D0 )))), TURNAS ) * DAS2R */
            L = ((T * (T * (T * (T * -2.447e-4 + .051635) + 31.8792) + 1717915923.2178) + 485868.249036) %
                Constants.TURNAS) * Constants.DAS2R;
            /*  Compute the fundamental argument LP. */
            /*    LP = MOD (       1287104.79305D0 +
                 .            T*( 129596581.0481D0 +
                 .            T*(       - 0.5532D0 +
                 .            T*(         0.000136D0 +
                 .            T*(       - 0.00001149D0 )))), TURNAS ) * DAS2R */
            LP = ((T * (T * (T * (T * -1.149e-5 + 1.36e-4) - .5532) + 129596581.0481) + 1287104.79305) %
                Constants.TURNAS) * Constants.DAS2R;
            /*  Compute the fundamental argument F. */
            /*    F  = MOD (       335779.526232D0 +
                 .                  T*( 1739527262.8478D0 +
                 .                  T*(       - 12.7512D0 +
                 .                  T*(       -  0.001037D0 +
                 .                  T*(          0.00000417D0 )))), TURNAS ) * DAS2R */
            F = ((T * (T * (T * (T * 4.17e-6 - .001037) - 12.7512) + 1739527262.8478) + 335779.526232) %
                Constants.TURNAS) * Constants.DAS2R;
            /*  Compute the fundamental argument D. */
            /*    D = MOD (        1072260.70369D0 +
                 .          T*( 1602961601.2090D0 +
                 .          T*(        - 6.3706D0 +
                 .          T*(          0.006593D0 +
                 .          T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R */
            D = ((T * (T * (T * (T * -3.169e-5 + .006593) - 6.3706) + 1602961601.209) + 1072260.70369) %
                Constants.TURNAS) * Constants.DAS2R;

            /*  Compute the fundamental argument OM. */
            /*    OM = MOD (       450160.398036D0 +
                 .             T*( - 6962890.5431D0 +
                 .             T*(         7.4722D0 +
                 .             T*(         0.007702D0 +
                 .             T*(       - 0.00005939D0 )))), TURNAS ) * DAS2R */
            OM = ((T * (T * (T * (T * -5.939e-5 + .007702) + 7.4722) - 6962890.5431) + 450160.398036) %
                Constants.TURNAS) * Constants.DAS2R;
        }


        private static class cnmtxD
        {
            /// <summary>
            /// tidal potential model for 71 diurnal and semidiurnal lines
            /// </summary>
            public static class tidalPotential
            {
                public static readonly double[] nj = new double[71] 
                {
                    2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
	                2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
	                2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
	                2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0
                };
                public static readonly double[] mj = new double[71] 
                {
                    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
                    1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,
	                1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,
	                2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0 
                };
                public static readonly double[] hs = new double[71]
                {
                    -1.94,-1.25,-6.64,-1.51,-8.02,-9.47,-50.2,
	                -1.8,-9.54,1.52,-49.45,-262.21,1.7,3.43,1.94,1.37,7.41,20.62,4.14,
	                3.94,-7.14,1.37,-122.03,1.02,2.89,-7.3,368.78,50.01,-1.08,2.93,
	                5.25,3.95,20.62,4.09,3.42,1.69,11.29,7.23,1.51,2.16,1.38,1.8,4.67,
	                16.01,19.32,1.3,-1.02,-4.51,120.99,1.13,22.98,1.06,-1.9,-2.18,
	                -23.58,631.92,1.92,-4.66,-17.86,4.47,1.97,17.2,294.0,-2.46,-1.02,
	                79.96,23.83,2.59,4.47,1.95,1.17
                };
                public static readonly double[] phase = new double[71]
                { 
                    9.0899831,8.8234208,12.1189598,1.44257,
	                4.738109,4.4715466,7.7670857,-2.9093042,.3862349,-3.1758666,
	                .1196725,3.4152116,12.8946194,5.5137686,6.4441883,-4.2322016,
	                -.9366625,8.5427453,11.8382843,1.1618945,5.9693878,-1.2032249,
	                2.0923141,-1.7847596,8.0679449,.8953321,4.1908712,7.4864102,
	                10.7819493,.3137975,6.2894282,7.2198478,-.161003,3.1345361,
	                2.8679737,-4.5128771,4.9665307,8.2620698,11.5576089,.6146566,
	                3.9101957,20.6617051,13.2808543,16.309831,8.9289802,5.0519065,
	                15.8350306,8.6624178,11.9579569,8.0808832,4.5771061,.7000324,
	                14.9869335,11.4831564,4.3105437,7.6060827,3.729009,10.6350594,
	                3.2542086,12.7336164,16.0291555,10.160259,6.2831853,2.4061116,
	                5.0862033,8.3817423,11.6772814,14.9728205,4.0298682,7.3254073,
	                9.1574019
                };
                public static readonly double[] freq = new double[71]
                {
                    5.1868805,5.38346657,5.38439079,5.41398343,
	                5.41490765,5.61149372,5.61241794,5.64201057,5.64293479,5.83859664,
	                5.83952086,5.84044508,5.84433381,5.87485066,6.03795537,6.06754801,
	                6.06847223,6.07236095,6.07328517,6.10287781,6.24878055,6.2650583,
	                6.26598252,6.28318449,6.28318613,6.29946388,6.3003881,6.30131232,
	                6.30223654,6.31759007,6.33479368,6.49789839,6.52841524,6.52933946,
	                6.72592553,6.75644239,6.76033111,6.76125533,6.76217955,6.98835826,
	                6.98928248,11.45675174,11.4872686,11.68477889,11.71529575,
	                11.73249771,11.89560406,11.91188181,11.91280603,11.930008,
	                11.94332289,11.96052486,12.11031632,12.12363121,12.13990896,
	                12.14083318,12.15803515,12.33834347,12.36886033,12.37274905,
	                12.37367327,12.54916865,12.56637061,12.58357258,12.59985198,
	                12.6007762,12.60170041,12.60262463,12.82880334,12.82972756,
	                13.06071921 
                };
            }

            /// <summary>
            /// Define the orthotide weight factors
            /// </summary>
            public static readonly double[,] sp = new double[6, 2]
            {
                {0.0298, 0.0200},
                {0.1408, 0.0905},
                {0.0805, 0.0638},
                {0.6002, 0.3476},
                {0.3025, 0.1645},
                {0.1517, 0.0923}
            };


        }

        public static double[] cnmtx(double dmjd) //       SUBROUTINE CNMTX ( DMJD,H )
        {
            /*  - - - - - - - - - -
             *   C N M T X
             *  - - - - - - - - - -

             *  This routine is part of the International Earth Rotation and
             *  Reference Systems Service (IERS) Conventions software collection.

             *  The purpose of the subroutine is to compute the time dependent part
             *  of second degree diurnal and semidiurnal tidal potential from
             *  the dominant spectral lines in the Cartwright-Tayler-Edden harmonic
             *  decomposition.

             *  In general, Class 1, 2, and 3 models represent physical effects that
             *  act on geodetic parameters while canonical models provide lower-level
             *  representations or basic computations that are used by Class 1, 2, or
             *  3 models.

             *  Status: Canonical model

             *     Class 1 models are those recommended to be used a priori in the
             *     reduction of raw space geodetic data in order to determine
             *     geodetic parameter estimates.
             *     Class 2 models are those that eliminate an observational
             *     singularity and are purely conventional in nature.
             *     Class 3 models are those that are not required as either Class
             *     1 or 2.
             *     Canonical models are accepted as is and cannot be classified as a
             *     Class 1, 2, or 3 model.

             *  Given:
             *     dmjd           d     Modified Julian Date

             *  Returned:
             *     h              d     vector of length 12 with partials of the
             *                          tidal variation with respect to the
             *                          orthoweights (Note 1)

             *  Notes:

             *  1) The diurnal and semidiurnal orthoweights fit to the 8 constituents
             *     are listed in Reference Ray et al.

             *  Test case:
             *     given input: dmjd = 54964.0D0

             *     expected output: h(1) = 15.35873641938967360D0
             *                      h(2) = 9.784941251812741214D0
             *                      h(3) = -5.520740128266865554D0
             *                      h(4) = 3.575314211234633888D0
             *                      h(5) = -13.93717453496387648D0
             *                      h(6) = -9.167400321705855504D0
             *                      h(7) = 5.532815475865292321D0
             *                      h(8) = 9.558741883500834646D0
             *                      h(9) = -10.22541212627272600D0
             *                      h(10)= 0.8367570529461261231D0
             *                      h(11)= 1.946355176475630611D0
             *                      h(12)= -13.55702062247304696D0

             *  References:

             *     Ray,R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
             *     "Diurnal and Semidiurnal Variations in the Earth's Rotation
             *     Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832

             *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
             *     IERS Technical Note No. 36, BKG (2010)

             *  Revisions:
             *  2008 November 07 B.E. Stetzler    Added header and copyright
             *  2008 November 21 B.E. Stetzler    Redefined variables as double
             *                                    precision and changed twopi
             *  2009 May      12 B.E. Stetzler    Added D0 to DATA structure values
             *  2009 May      19 B.E. Stetzler    Provided validated test case
             *  2009 June     08 B.E. Stetzler    Redefined pinm and alpha, used
             *                                    double precision exclusively, and
             *                                    updated validated test case values
             *  2009 September 03 B.E.Stetzler    Capitalized all variables for Fortran
             *                                    77 compatibility
             *  2010 March     17 B.E.Stetzler    Aligned table array values in DATA
             *                                    structures
            /* ----------------------------------------------------------------------- */

            const double d1960 = 37076.5;
            const double dt = 2.0;
            const int nmax = 2;

            double[, ,] anm = new double[2, 4, 3];
            double[, ,] bnm = new double[2, 4, 3];

            const int NLINES = 71;

            /* Compute the time dependent potential matrix */
            for (int k = -1; k <= 1; k++)                           //       DO 50 K=-1,1
            {
                double dt60 = dmjd - k * dt - d1960;                //          DT60 = (DMJD - K*DT) - D1960

                for (int j = 0; j < NLINES; j++)                    //          DO 100 J=1,NLINES
                {
                    int n = (int)cnmtxD.tidalPotential.nj[j];                             //             N = NJ(J)
                    int m = (int)cnmtxD.tidalPotential.mj[j];                             //             M = MJ(J)

                    double pinm = (double)((n + m) % 2) * Constants.TWOPI / 4.0;     //             PINM = DFLOAT(MOD(N+M,2))*TWOPI/4.0D0
                    /* alpha = phase(j) + freq(j)*dt60 - pinm */
                    double alpha = ((cnmtxD.tidalPotential.phase[j] - pinm) % Constants.TWOPI) +
                                   ((cnmtxD.tidalPotential.freq[j] * dt60) % Constants.TWOPI);        //            ALPHA = DMOD(PHASE(J) - PINM,TWOPI) + DMOD(FREQ(J)*DT60,TWOPI)
                    anm[n - 2, m, k + 1] += cnmtxD.tidalPotential.hs[j] * Math.Cos(alpha);    //             ANM(N,M,K) = ANM(N,M,K) + HS(J)*DCOS(ALPHA)
                    bnm[n - 2, m, k + 1] -= cnmtxD.tidalPotential.hs[j] * Math.Sin(alpha);    //             BNM(N,M,K) = BNM(N,M,K) - HS(J)*DSIN(ALPHA)
                }                                                   // 100      CONTINUE
            }                                                       // 50    CONTINUE

            /* orthogonalize the response terms */
            double[,] p = new double[3, 2];
            double[,] q = new double[3, 2];

            for (int m = 1; m < 3; m++)                             //       DO 150 M = 1,2
            {
                double ap = anm[0, m, 2] + anm[0, m, 0];                   //         AP = ANM(2,M,1) + ANM(2,M,-1)
                double am = anm[0, m, 2] - anm[0, m, 0];                   //         AM = ANM(2,M,1) - ANM(2,M,-1)
                double bp = bnm[0, m, 2] + bnm[0, m, 0];                   //         BP = BNM(2,M,1) + BNM(2,M,-1)
                double bm = bnm[0, m, 2] - bnm[0, m, 0];                   //         BM = BNM(2,M,1) - BNM(2,M,-1)
                p[0, m - 1] = cnmtxD.sp[0, m - 1] * anm[0, m, 1];  //         P(0,M) = SP(1,M)*ANM(2,M,0)
                p[1, m - 1] = cnmtxD.sp[1, m - 1] * anm[0, m, 1] -
                              cnmtxD.sp[2, m - 1] * ap;            //         P(1,M) = SP(2,M)*ANM(2,M,0) - SP(3,M)*AP
                p[2, m - 1] = cnmtxD.sp[3, m - 1] * anm[0, m, 1] -
                              cnmtxD.sp[4, m - 1] * ap +
                              cnmtxD.sp[5, m - 1] * bm;            //         P(2,M) = SP(4,M)*ANM(2,M,0) - SP(5,M)*AP + SP(6,M)*BM

                q[0, m - 1] = cnmtxD.sp[0, m - 1] * bnm[0, m, 1];  //         Q(0,M) = SP(1,M)*BNM(2,M,0)
                q[1, m - 1] = cnmtxD.sp[1, m - 1] * bnm[0, m, 1] -
                              cnmtxD.sp[2, m - 1] * bp;            //         Q(1,M) = SP(2,M)*BNM(2,M,0) - SP(3,M)*BP
                q[2, m - 1] = cnmtxD.sp[3, m - 1] * bnm[0, m, 1] -
                              cnmtxD.sp[4, m - 1] * bp -
                              cnmtxD.sp[5, m - 1] * am;            //         Q(2,M) = SP(4,M)*BNM(2,M,0) - SP(5,M)*BP - SP(6,M)*AM
                anm[0, m, 0] = p[0, m - 1];                             //         ANM(2,M,-1) = P(0,M)
                anm[0, m, 1] = p[1, m - 1];                             //         ANM(2,M,0) = P(1,M)
                anm[0, m, 2] = p[2, m - 1];                             //         ANM(2,M,1) = P(2,M)
                bnm[0, m, 0] = q[0, m - 1];                             //         BNM(2,M,-1) = Q(0,M)
                bnm[0, m, 1] = q[1, m - 1];                             //         BNM(2,M,0) = Q(1,M)
                bnm[0, m, 2] = q[2, m - 1];                             //         BNM(2,M,1) = Q(2,M)
            }                                                       // 150   CONTINUE

            /* fill partials vector */
            double[] h = new double[12];
            int jj = 0;                                             //       J = 1
            for (int n = 2; n <= nmax; n++)                         //       DO 200 N=2,NMAX
            {
                for (int m = 1; m <= n; m++)                        //          DO 250 M = 1,N
                {
                    for (int k = -1; k <= 1; k++)                   //             DO 300 K = -1,1
                    {
                        h[jj] = anm[n - 2, m, k + 1];                       //                H(J)  = ANM(N,M,K)
                        h[jj + 1] = bnm[n - 2, m, k + 1];                   //                H(J+1)= BNM(N,M,K)
                        jj += 2;                                    //                J = J + 2
                    }                                               // 300         CONTINUE
                }                                                   // 250      CONTINUE 
            }                                                       // 200   CONTINUE
            return h;                                               //       RETURN
        }
    }
}
