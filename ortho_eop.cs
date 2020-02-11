using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IERS_Routines
{
    public class OceanTides
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
        /// Struct returning the calculated values, X-pol, Y-Pol and UT1 variations according to the ocean tide effect.
        /// </summary>
        public struct PM
        {
            public double dX;
            public double dY;
            public double dUT1;
        }

        #region public parameter
        /// <summary>
        /// Get/Set the Unit for x pol variation (has to be an angle unit, according to the unit class)
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_X
        {
            get
            {
                return _OutputUnit_X;
            }
            set
            {
                _OutputUnit_X = value;
                dXc = Units.Units.getConversion(BaseUnits.X, _OutputUnit_X).Factor;
            }
        }

        /// <summary>
        /// Internal unit for the x pol variations to calculate the conversion factor dX.
        /// </summary>
        internal Units.Units.UnitNamesEnum _OutputUnit_X = BaseUnits.X;

        /// <summary>
        /// Internal conversion factor for the x pol variations.
        /// </summary>
        internal double dXc = 1;

        /// <summary>
        /// Get/Set the Unit for y pol variation (has to be an angle unit, according to the unit class)
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_Y
        {
            get
            {
                return _OutputUnit_Y;
            }
            set
            {
                _OutputUnit_Y = value;
                dYc = Units.Units.getConversion(BaseUnits.Y, _OutputUnit_Y).Factor;
            }
        }

        /// <summary>
        /// Internal unit for the y pol variations to calculate the conversion factor dY.
        /// </summary>
        internal Units.Units.UnitNamesEnum _OutputUnit_Y = BaseUnits.Y;

        /// <summary>
        /// Internal conversion factor for the y pol variations.
        /// </summary>
        internal double dYc = 1;

        /// <summary>
        /// Get/Set the Unit for UT1 variation (has to be a time unit, according to the unit class)
        /// </summary>
        public Units.Units.UnitNamesEnum OutputUnit_UT1
        {
            get
            {
                return _OutputUnit_UT1;
            }
            set
            {
                _OutputUnit_UT1 = value;
                dUT1c = Units.Units.getConversion(BaseUnits.UT1, _OutputUnit_UT1).Factor;
            }
        }

        /// <summary>
        /// Internal unit for the UT1 variations to calculate the conversion factor dUT1.
        /// </summary>
        internal Units.Units.UnitNamesEnum _OutputUnit_UT1 = BaseUnits.UT1;

        /// <summary>
        /// Internal conversion factor for the UT1 variations.
        /// </summary>
        internal double dUT1c = 1;
        #endregion

        #region base units
        public static class BaseUnits
        {
            public readonly static Units.Units.UnitNamesEnum X = Units.Units.UnitNamesEnum.MicroArcSecond;
            public readonly static Units.Units.UnitNamesEnum Y = Units.Units.UnitNamesEnum.MicroArcSecond;
            public readonly static Units.Units.UnitNamesEnum UT1 = Units.Units.UnitNamesEnum.MicroSecond;
        }
        #endregion

        public OceanTides()
        { }

        /// <summary>
        /// diurnal and semidiurnal orthoweights fit to the 8 constituents
        /// listed in IERS Technical Note 21, July 1996 which are from the
        /// paper "Diurnal and Semidiurnal Variations in the Earth's Rotation
        /// Rate Induced by Ocean Tides"  by Ray,R.D., Steinberg,D.J.,
        /// Chao,B.F., and Cartwright, D.E., Science, 264, pp. 830-832.
        /// </summary>
        private static readonly double[,] orthow = new double[12, 3]  
            {
                {  -6.77832, 14.86283, -1.76335},   
                { -14.86323, -6.77846,  1.03364},   
                {   0.47884,  1.45234, -0.27553},   
                {  -1.45303,  0.47888,  0.34569},   
                {   0.16406, -0.42056, -0.12343},   
                {   0.42030,  0.16469, -0.10146},   
                {   0.09398, 15.30276, -0.47119},   
                {  25.73054, -4.30615,  1.28997},   
                {  -4.77974,  0.07564, -0.19336},   
                {   0.28080,  2.28321,  0.02724},   
                {   1.94539, -0.45717,  0.08955},   
                {  -0.73089, -1.62010,  0.04726}  
            };

        public PM ortho_eop(double time)        //       SUBROUTINE ORTHO_EOP ( TIME, EOP )
        {
            /*  - - - - - - - - - -
             *   O R T H O _ E O P
             *  - - - - - - - - - -

             *  This routine is part of the International Earth Rotation and
             *  Reference Systems Service (IERS) Conventions software collection.

             *  The purpose of the subroutine is to compute the diurnal and semi-
             *  diurnal variations in Earth Orientation Parameters (x,y, UT1) from
             *  ocean tides.

             *  In general, Class 1, 2, and 3 models represent physical effects that
             *  act on geodetic parameters while canonical models provide lower-level
             *  representations or basic computations that are used by Class 1, 2, or
             *  3 models.

             *  Status: Class 1 model

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
             *     TIME           d     Modified Julian Date

             *  Returned:
             *     EOP            d     delta_x, in microarcseconds
             *                          delta_y, in microarcseconds
             *                          delta_UT1, in microseconds

             *  Notes:

             *  1) The diurnal and semidiurnal orthoweights fit to the 8 constituents
             *     are listed in Reference 1.

             *  Called:
             *     CNMTX                Compute time dependent part of second degree
             *                          diurnal and semidiurnal tidal potential from
             *                          the dominant spectral lines in the Cartwright-
             *                          Tayler-Edden harmonic decomposition
             *  Test case:
             *     given input: MJD = 47100D0

             *     expected output: delta_x = -162.8386373279636530D0 microarcseconds
             *                      delta_y = 117.7907525842668974D0 microarcseconds
             *                      delta_UT1 = -23.39092370609808214D0 microseconds

             *  References:

             *     Ray, R. D., Steinberg, D. J., Chao, B. F., and Cartwright, D. E.,
             *     "Diurnal and Semidiurnal Variations in the Earth's Rotation
             *     Rate Induced by Ocean Tides", 1994, Science, 264, pp. 830-832

             *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
             *     IERS Technical Note No. 36, BKG (2010)

             *  Revisions:
             *  1997 February    R. Eanes         Original code
             *  2008 November 07 B.E. Stetzler    Added header and copyright
             *  2008 November 07 B.E. Stetzler    Provided test case
             *  2009 May      12 B.E. Stetzler    Replaced ENDDO statements with
             *                                    CONTINUE statements
             *  2009 June     09 B.E. Stetzler    Updated validated test case values
             *                                    based on changes to CNMTX.F
             *  2010 March    19 B.E. Stetzler    Capitalized variables for FORTRAN
             *                                    77 backwards compatibility
             * ----------------------------------------------------------------------- */

            /* Compute the partials of the tidal variations to the orthoweights */
            double[] h = CommonFunctions.cnmtx(time);   //       CALL CNMTX (TIME, H)
            
            /* Compute eop changes */
            double[] eop = new double[3];
            for (int k = 0; k < 3; k++)                 //       DO 20 K=1,3
            {
                eop[k] = 0.0;                           //          EOP(K) = 0D0
                for (int j = 0; j < 12; ++j)            //          DO 40 J=1,12
                {
                    eop[k] += (h[j] * orthow[j, k]);  //             EOP(K) = EOP(K) + H(J)*ORTHOW(J,K)
                }                                       // 40       CONTINUE
            }                                           // 20    CONTINUE

            return new PM
            { 
                dX = eop[0] * dXc,
                dY = eop[1] * dYc,
                dUT1 = eop[2] * dUT1c
            };       //       RETURN
        }




        public PM ortho_eop(double mjd, double stationLongitude)
        {
            double stlong = Math.PI / 180d * stationLongitude;

            double sa = Math.Sin(stlong);
            double ca = Math.Cos(stlong);

            PM pm = ortho_eop(mjd);

            return new PM
            {
                dX = pm.dX * ca - pm.dY * sa,
                dY = -pm.dX * sa + pm.dY * ca,
                dUT1 = pm.dUT1,
            };
        }
    }
}
