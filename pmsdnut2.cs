using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IERS_Routines
{
    public class NutationModel
    {
        /*  Copyright (C) 2008
         *  IERS Conventions Center

         *  ==================================
         *  IERS Conventions Software License
         *  ==================================

         *  NOTICE TO USER:

         *  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
         *  WHICH APPLY TO ITS USE.

         *  1. The Software is provided by the IERS Conventions Center ("the
         *     Center").

         *  2. Permission is granted to anyone to use the Software for any
         *     purpose, including commercial applications, free of charge,
         *     subject to the conditions and restrictions listed below.

         *  3. You (the user) may adapt the Software and its algorithms for your
         *     own purposes and you may distribute the resulting "derived work"
         *     to others, provided that the derived work complies with the
         *     following requirements:

         *     a) Your work shall be clearly identified so that it cannot be
         *        mistaken for IERS Conventions software and that it has been
         *        neither distributed by nor endorsed by the Center.

         *     b) Your work (including source code) must contain descriptions of
         *        how the derived work is based upon and/or differs from the
         *        original Software.

         *     c) The name(s) of all modified routine(s) that you distribute
         *        shall be changed.

         *     d) The origin of the IERS Conventions components of your derived
         *        work must not be misrepresented; you must not claim that you
         *        wrote the original Software.

         *     e) The source code must be included for all routine(s) that you
         *        distribute.  This notice must be reproduced intact in any
         *        source distribution.

         *  4. In any published work produced by the user and which includes
         *     results achieved by using the Software, you shall acknowledge
         *     that the Software was used in obtaining those results.

         *  5. The Software is provided to the user "as is" and the Center makes
         *     no warranty as to its use or performance.   The Center does not
         *     and cannot warrant the performance or results which the user may
         *     obtain by using the Software.  The Center makes no warranties,
         *     express or implied, as to non-infringement of third party rights,
         *     merchantability, or fitness for any particular purpose.  In no
         *     event will the Center be liable to the user for any consequential,
         *     incidental, or special damages, including any lost profits or lost
         *     savings, even if a Center representative has been advised of such
         *     damages, or for any claim by any third party.

         *  Correspondence concerning IERS Conventions software should be
         *  addressed as follows:

         *                     Gerard Petit
         *     Internet email: gpetit[at]bipm.org
         *     Postal address: IERS Conventions Center
         *                     Time, frequency and gravimetry section, BIPM
         *                     Pavillon de Breteuil
         *                     92312 Sevres  FRANCE
         *     or
         *                     Brian Luzum
         *     Internet email: brian.luzum[at]usno.navy.mil
         *     Postal address: IERS Conventions Center
         *                     Earth Orientation Department
         *                     3450 Massachusetts Ave, NW
         *                     Washington, DC 20392
         * -----------------------------------------------------------------------*/

        /// <summary>
        /// Contructor for calculation the nutation model, adjustable units
        /// </summary>
        public NutationModel()
        { }

        /// <summary>
        /// Struct containing calculated pole coordinates
        /// </summary>
        public struct PM
        {
            /// <summary>
            /// Polar coordinate X (oriented along Greenwich meridian)
            /// </summary>
            public double dX;

            /// <summary>
            /// Polar coordinate Y (oriented along 90 degree west)
            /// </summary>
            public double dY;
        }

        /// <summary>
        /// Gets/Sets the output unit for calculation.
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit
        {
            get
            {
                return _Unit;
            }
            set
            {
                _Unit = value;
                dConversion = Units.Units.getConversion(BaseUnit, _Unit).Factor;
            }
        }

        /// <summary>
        /// Internal parmater for selection and setting the output unit.
        /// </summary>
        internal Units.Units.UnitNamesEnum _Unit = BaseUnit;

        /// <summary>
        /// Internal parmater for setting conversion factor to the output unit.
        /// </summary>
        internal double dConversion = 1d;

        /// <summary>
        /// Base unit, where the vales nativly calculated. 
        /// </summary>
        public static readonly Units.Units.UnitNamesEnum BaseUnit = Units.Units.UnitNamesEnum.MicroArcSecond;

        public static class data
        {
            public static readonly int[,] iarg = new int[25, 6]
            {
                { 0,  0, 0,  0,  0, -1 },   // Coefficients of the long periodic terms in polar motion
                { 0, -1, 0,  1,  0,  2 },   // Source: IERS Conventions (2010), Table 5.1a   
                { 0, -1, 0,  1,  0,  1 },
                { 0, -1, 0,  1,  0,  0 },
                { 0,  1, 1, -1,  0,  0 },
                { 0,  1, 1, -1,  0, -1 },
                { 0,  0, 0,  1, -1,  1 },
                { 0,  1, 0,  1, -2,  1 },
                { 0,  0, 0,  1,  0,  2 },
                { 0,  0, 0,  1,  0,  1 },
                { 0,  0, 0,  1,  0,  0 },
                { 0, -1, 0,  1,  2,  1 },
                { 0,  1, 0,  1,  0,  1 },
                { 0,  0, 0,  3,  0,  3 },
                { 0,  0, 0,  3,  0,  2 },
                { 1, -1, 0, -2,  0, -1 },   // Coefficients of the quasi diurnal terms in polar motion
                { 1, -1, 0, -2,  0, -2 },   // Source: IERS Conventions (2010), Table 5.1a
                { 1,  1, 0, -2, -2, -2 },
                { 1,  0, 0, -2,  0, -1 },
                { 1,  0, 0, -2,  0, -2 },
                { 1, -1, 0,  0,  0,  0 },
                { 1,  0, 0, -2,  2, -2 },
                { 1,  0, 0,  0,  0,  0 },
                { 1,  0, 0,  0,  0, -1 },
                { 1,  1, 0,  0,  0,  0 },
            };
            public static readonly double[] xs = new double[25] { 0.0, 1.5, -28.5, -4.7, -0.7, 1.0, 1.2, 1.3, -0.1, 0.9, 0.1, 0.0, -0.1, -0.1, -0.1, -0.4, -2.3, -0.4, -2.1, -11.4, 0.8, -4.8, 14.3, 1.9, 0.8 };
            public static readonly double[] xc = new double[25] { 0.6, 0.0, -0.2, -0.1, 0.2, 0.3, 0.2, 0.4, -0.2, 4.0, 0.6, 0.1, 0.3, 0.1, 0.1, 0.3, 1.3, 0.3, 1.2, 6.5, -0.5, 2.7, -8.2, -1.1, -0.4 };
            public static readonly double[] ys = new double[25] { -0.1, -0.2, 3.4, 0.6, -0.2, -0.3, -0.2, -0.2, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, -0.3, -1.3, -0.3, -1.2, -6.5, 0.5, -2.7, 8.2, 1.1, 0.4 };
            public static readonly double[] yc = new double[25] { -0.1, 0.1, -3.9, -0.9, -0.7, 1.0, 1.4, 2.9, -1.7, 32.4, 5.1, 0.6, 2.7, 0.9, 0.6, -0.4, -2.3, -0.4, -2.1, -11.4, 0.8, -4.8, 14.3, 1.9, 0.8 };

            /* Rate of secular polar motion, in microarcseconds per year */
            /* Source: IERS Conventions (2010), Table 5.1a */
            public static readonly double xrate = -3.8;
            public static readonly double yrate = -4.3;
        }

        /// <summary>
        /// Evaluates the model of polar motion for a nonrigid Earth due to tidal gravitation
        /// </summary>
        /// <param name="rmjd">Epoch in modified julian day for calculation of polar motion </param>
        /// <param name="IBAND">Parameter defining the range of periods for the terms which are included in computations; if equal to 1 only the quasi diurnal terms are computed, otherwise the full model</param>
        /// <returns></returns>
        public PM pmsdnut2(double rmjd, bool IBAND = true)     // SUBROUTINE PMSDNUT2 (RMJD, PM)
        {
           /*  - - - - - - - - - - -
            *   P M S D N U T 2
            *  - - - - - - - - - - -

            *  This routine is part of the International Earth Rotation and
            *  Reference Systems Service (IERS) Conventions software collection.

            *  This subroutine evaluates the model of polar motion for
            *  a nonrigid Earth due to tidal gravitation. This polar motion
            *  is equivalent to the so-called "subdiurnal nutation." The model
            *  is a sum of a first order polynomial and 25 trigonometric terms
            *  (15 long periodic and 10 quasi diurnal) with coefficients given
            *  in Table 5.1a of the IERS Conventions (2010).

            *     :------------------------------------------:
            *     :                                          :
            *     :                 IMPORTANT                :
            *     :                                          :
            *     : In the present version this subroutine   :
            *     : neglects the linear trend and the long   :
            *     : periodic terms of the expansion, for the :
            *     : reasons explained in Section 5.x.x of    :
            *     : the IERS Conventions (2010), last para-  :
            *     : graph before Table 5.1. If the full      :
            *     : expansion is needed, set the parameter   :
            *     : iband to 0 instead of 1, that is replace :
            *     : the statement                            :
            *     :     PARAMETER ( iband = 1 )              :
            *     : to  PARAMETER ( iband = 0 )              :
            *     :                                          :
            *     :__________________________________________:

            *  In general, Class 1, 2, and 3 models represent physical effects that
            *  act on geodetic parameters while canonical models provide lower-level
            *  representations or basic computations that are used by Class 1, 2, or
            *  3 models.

            *  Status:  Class 3 model

            *     Class 1 models are those recommended to be used a priori in the
            *     reduction of raw space geodetic data in order to determine
            *     geodetic parameter estimates.
            *     Class 2 models are those that eliminate an observational
            *     singularity and are purely conventional in nature.
            *     Class 3 models are those that are not required as either Class
            *     1 or 2.
            *     Canonical models are accepted as is and cannot be classified as
            *     a Class 1, 2, or 3 model.

            *  Given:
            *     rmjd        d      Time expressed as modified Julian date

            *  Returned:
            *     pm          d(2)      Vector of length 2 (Note 1)

            *  Notes:

            *  1) The polar motion coordinates (dx, dy) are expressed in
            *     microarcseconds.

            *  Called:
            *     FUNDARG             Compute the angular fundamental arguments

            *  Test case:
            *     given input: rmjd = 54335D0 ( August 23, 2007 )

            *     expected output: (dx) pm(1)  = 24.65518398386097942D0 microarcseconds
            *                      (dy) pm(2) = -14.11070254891893327D0 microarcseconds

            *  References:

            *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
            *     IERS Technical Note No. 36, BKG (2010)

            *  Revisions:
            *  2005 March       A.Brzezinski   Original code
            *  2008 November 26 B.E.Stetzler   Initial changes to code
            *  2008 December 01 B.E.Stetzler   Provided test case
            *  2009 August   18 B.E.Stetzler   Capitalized all variables for FORTRAN
            *                                  77 compatibility
            *  2010 May      14 B.E.Stetzler   Replaced call to PMARGS to FUNDARG
            *                                  for universal fundamental argument
            *                                  subroutine
            *  2010 May      17 B.E.Stetzler   Validated test case using internally
            *                                  computed GMST and call to FUNDARG
            *                                  matched previous external call to
            *                                  PMARGS for all six parameters
            *  2010 June     23 B.E.Stetzler   Modified coefficients of the long
            *                                  and short period terms in polar
            *                                  motion and secular polar motion
            *                                  rate to coincide with Table 5.1a
            *         ----------------------------
            *           D E F I N I T I O N S
            *         ----------------------------
            *  iband  - parameter defining the range of periods for the terms which
            *           are included in computations; if equal to 1 only the quasi
            *           diurnal terms are computed, otherwise the full model
            *  iarg   - array defining for each of the 25 trigonometric terms a set
            *           of 6 integer multipliers of the fundamental angular arguments
            *  arg    - vector of the following 6 fundamental arguments used to
            *           compute the angular argument of the trigonometric functions
            *           arg(1:6) = [ GMST+pi, el, elp, f, d, om ]; this vector is
            *           evaluated by the subroutine FUNDARG which is called as an
            *           external subroutine.  Originally evaluated by the subroutine
            *           PMARGS.
            *  period - array of periods of the trigonometric terms of expansion, in
            *           mean solar days; only for a check - not used in computations
            *  xs, xc - sine and cosine coefficients of the x coordinate of the pole,
            *           in microarcseconds
            *  ys, yc - sine and cosine coefficients of the y coordinate of the pole,
            *           in microarcseconds
            *  angle  - angular argument of the trigonometric functions
            *           angle = Sum(i=1:6) iarg(i,j)*arg(i), for j=1,25 */


            /* Compute the periodical part of the model */

            /* Coordinates of the pole are set to zero first */
            PM pm = new PM { dX = 0.0, dY = 0.0 };              //       PM(1) = 0D0 & PM(2) = 0D0

            /* Evaluate the vector of the fundamental arguments 
             * arg(1:6) = [ GMST+pi, el, elp, f, d, om ] at t = rmjd */

            /*  Convert the input epoch to Julian centuries of TDB since J2000 */
            double T = (rmjd - Constants.RMJD0) / 36525.0;                //       T = (RMJD-RMJD0)/36525D0 

            /*  Compute GMST + pi */
            double GMST = (T * (T * (T * -6.2e-6 + .093104) +
                3164400184.8128662) + 67310.54841) % 86400.0;  /*       GMST = MOD (   67310.54841D0 +
                                                                       .               T*( (8640184.812866D0 + 3155760000D0) +
                                                                       .               T*( 0.093104D0 +
                                                                       .               T*( -0.0000062 ))), 86400D0 ) */

            double L = 0.0, F = 0.0, D = 0.0, LP = 0.0, OM = 0.0;

            CommonFunctions.fundarg(T, ref L, ref LP, ref F, ref D, ref OM);    //       CALL FUNDARG ( T, L, LP, F, D, OM )

            double[] arg = new double[6];
            arg[0] = GMST / Constants.RAD2SEC + Constants.PI;                       //       ARG(1) = GMST / RAD2SEC + PI
            arg[0] = arg[0] % Constants.TWOPI;                            //       ARG(1) = DMOD( ARG(1), TWOPI )
            arg[1] = L;                                         //       ARG(2) = L
            arg[2] = LP;                                        //       ARG(3) = LP
            arg[3] = F;                                         //       ARG(4) = F
            arg[4] = D;                                         //       ARG(5) = D
            arg[5] = OM;                                        //       ARG(6) = OM 

            int jstart;
            if (IBAND)                                     //       IF (IBAND.EQ.1) THEN
                jstart = 15;                                    //         JSTART = 16
            else                                                //       ELSE
                jstart = 0;                                     //         JSTART = 1

            for (int j = jstart; j < 25; ++j)                   //       DO 20 J=JSTART,25
            {
                /* For the j-th term of the trigonometric expansion, compute the angular
                 * argument angle of sine and cosine functions as a linear integer
                 * combination of the 6 fundamental arguments */
                double angle = 0.0;                             //         ANGLE = 0D0
                for (int i = 0; i < 6; i++)                     //         DO 10 I=1,6
                {
                    angle += data.iarg[j, i] * arg[i];               //           ANGLE = ANGLE + IARG(I,J) * ARG(I)
                }
                angle = angle % Constants.TWOPI;                          //         ANGLE = DMOD( ANGLE, TWOPI )

                /* Compute contribution from the j-th term to the polar motion coordinates */
                double sa = Math.Sin(angle);
                double ca = Math.Cos(angle);
                pm.dX += data.xs[j] * sa + data.xc[j] * ca;               //         PM(1) = PM(1) + XS(J)*DSIN(ANGLE) + XC(J)*DCOS(ANGLE)
                pm.dY += data.ys[j] * sa + data.yc[j] * ca;               //         PM(2) = PM(2) + YS(J)*DSIN(ANGLE) + YC(J)*DCOS(ANGLE)
            }                                                   //    20 CONTINUE

            if (IBAND)                                     //       IF (IBAND.EQ.1) RETURN
                return new PM { dX = pm.dX * dConversion, dY = pm.dY * dConversion };

            /* Add the secular term of the model */
            pm.dX += data.xrate * (rmjd - Constants.RMJD0) / 365.25;           //       PM(1) = PM(1) + XRATE * (RMJD-RMJD0) / 365.25D0
            pm.dY += data.yrate * (rmjd - Constants.RMJD0) / 365.25;           //       PM(2) = PM(2) + YRATE * (RMJD-RMJD0) / 365.25D0

            return new PM { dX = pm.dX * dConversion, dY = pm.dY * dConversion };                                          //       RETURN
        }
    }
}
