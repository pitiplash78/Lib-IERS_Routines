using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IERS_Routines
{
    public class ZonalTides
    {
        /// <summary>
        /// Gets/Sets the component to be calculated.
        /// </summary>
        public Components Component { get; private set; }

        /// <summary>
        /// Enumerator class for selection of the component to be calculated.
        /// </summary>
        public enum Components
        {
            /// <summary>
            /// correction on UT1
            /// </summary>
            UT1,

            /// <summary>
            /// correction on LOD
            /// </summary>
            LOD,

            /// <summary>
            /// correction on Omega
            /// </summary>
            Omega,

            /// <summary>
            /// correction on dC/C (C is the axial moment of inertia of the whole Earth (8.037D37 kg.m2.s-1))
            /// </summary>
            dCC,
        }

        #region Base Units
        /// <summary>
        /// Base unit of related to the selected component to be calculated.
        /// </summary>
        public Units.Units.UnitNamesEnum BaseUnit;

        /// <summary>
        /// Defines the native calculation units, of according component.
        /// </summary>
        public static class BaseUnits
        {
            public static readonly Units.Units.UnitNamesEnum UT1 = Units.Units.UnitNamesEnum.Second;
            public static readonly Units.Units.UnitNamesEnum LOD = Units.Units.UnitNamesEnum.SecondPerDay;
            public static readonly Units.Units.UnitNamesEnum Omega = Units.Units.UnitNamesEnum.RadianPerSecond;
            public static readonly Units.Units.UnitNamesEnum dCC = Units.Units.UnitNamesEnum.Radian;
        }
        #endregion

        /// <summary>
        /// Gets/Sets the output unit and according internal parameter.
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
                dConversion = Units.Units.getConversion(this.BaseUnit, _Unit).Factor;
            }

        }

        /// <summary>
        /// Internal parameter of the unit, in which has to be calculated the according component.
        /// </summary>
        internal Units.Units.UnitNamesEnum _Unit = Units.Units.UnitNamesEnum.Unknown;

        /// <summary>
        /// Internal parameter for conversion the unit, in order to calculate the component.
        /// </summary>
        internal double dConversion = 1d;

        #region Data
        /// <summary>
        /// Tables of Luni-Solar arguments, multiples and coefficients 
        /// </summary>
        private static class Table
        {
            //  Luni-Solar argument multipliers

            public static readonly int NZONT = 62;

            /// <summary>
            /// Coefficients for the fundamental arguments
            /// </summary>
            public static readonly int[,] nfund = new int[62, 5]
            {
                // l,  l', F,  D,  OMEGA
                 { 1,  0,  2,  2,  2 }, 
                 { 2,  0,  2,  0,  1 }, 
                 { 2,  0,  2,  0,  2 }, 
                 { 0,  0,  2,  2,  1 }, 
                 { 0,  0,  2,  2,  2 }, 
                 { 1,  0,  2,  0,  0 }, 
                 { 1,  0,  2,  0,  1 }, 
                 { 1,  0,  2,  0,  2 }, 
                 { 3,  0,  0,  0,  0 }, 
                 {-1,  0,  2,  2,  1 }, 
                 {-1,  0,  2,  2,  2 }, 
                 { 1,  0,  0,  2,  0 }, 
                 { 2,  0,  2, -2,  2 }, 
                 { 0,  1,  2,  0,  2 }, 
                 { 0,  0,  2,  0,  0 }, 
                 { 0,  0,  2,  0,  1 }, 
                 { 0,  0,  2,  0,  2 }, 
                 { 2,  0,  0,  0, -1 }, 
                 { 2,  0,  0,  0,  0 }, 
                 { 2,  0,  0,  0,  1 }, 
                 { 0, -1,  2,  0,  2 }, 
                 { 0,  0,  0,  2, -1 }, 
                 { 0,  0,  0,  2,  0 }, 
                 { 0,  0,  0,  2,  1 }, 
                 { 0, -1,  0,  2,  0 }, 
                 { 1,  0,  2, -2,  1 }, 
                 { 1,  0,  2, -2,  2 }, 
                 { 1,  1,  0,  0,  0 }, 
                 {-1,  0,  2,  0,  0 }, 
                 {-1,  0,  2,  0,  1 }, 
                 {-1,  0,  2,  0,  2 }, 
                 { 1,  0,  0,  0, -1 }, 
                 { 1,  0,  0,  0,  0 }, 
                 { 1,  0,  0,  0,  1 }, 
                 { 0,  0,  0,  1,  0 }, 
                 { 1, -1,  0,  0,  0 }, 
                 {-1,  0,  0,  2, -1 }, 
                 {-1,  0,  0,  2,  0 }, 
                 {-1,  0,  0,  2,  1 }, 
                 { 1,  0, -2,  2, -1 }, 
                 {-1, -1,  0,  2,  0 }, 
                 { 0,  2,  2, -2,  2 }, 
                 { 0,  1,  2, -2,  1 }, 
                 { 0,  1,  2, -2,  2 }, 
                 { 0,  0,  2, -2,  0 }, 
                 { 0,  0,  2, -2,  1 }, 
                 { 0,  0,  2, -2,  2 }, 
                 { 0,  2,  0,  0,  0 }, 
                 { 2,  0,  0, -2, -1 }, 
                 { 2,  0,  0, -2,  0 }, 
                 { 2,  0,  0, -2,  1 }, 
                 { 0, -1,  2, -2,  1 }, 
                 { 0,  1,  0,  0, -1 }, 
                 { 0, -1,  2, -2,  2 }, 
                 { 0,  1,  0,  0,  0 }, 
                 { 0,  1,  0,  0,  1 }, 
                 { 1,  0,  0, -1,  0 }, 
                 { 2,  0, -2,  0,  0 }, 
                 {-2,  0,  2,  0,  1 }, 
                 {-1,  1,  0,  1,  0 },
                 { 0,  0,  0,  0,  2 }, 
                 { 0,  0,  0,  0,  1 }
            };

            /// <summary>
            /// Zonal tide term coefficients
            /// </summary>
            public static readonly double[,] tide = new double[62, 6] 
            {   
                /*    Multiple of
                 *    DUT             DLOD              DOMEGA
                 *    sin     cos      cos      sin       cos      sin */
                {    -0.0235, 0.0000,  0.2617,  0.0000,  -0.2209,  0.0000},
                {    -0.0404, 0.0000,  0.3706,  0.0000,  -0.3128,  0.0000},
                {    -0.0987, 0.0000,  0.9041,  0.0000,  -0.7630,  0.0000},
                {    -0.0508, 0.0000,  0.4499,  0.0000,  -0.3797,  0.0000},
                {    -0.1231, 0.0000,  1.0904,  0.0000,  -0.9203,  0.0000},
                {    -0.0385, 0.0000,  0.2659,  0.0000,  -0.2244,  0.0000},
                {    -0.4108, 0.0000,  2.8298,  0.0000,  -2.3884,  0.0000},
                {    -0.9926, 0.0000,  6.8291,  0.0000,  -5.7637,  0.0000},
                {    -0.0179, 0.0000,  0.1222,  0.0000,  -0.1031,  0.0000},
                {    -0.0818, 0.0000,  0.5384,  0.0000,  -0.4544,  0.0000},
                {    -0.1974, 0.0000,  1.2978,  0.0000,  -1.0953,  0.0000},
                {    -0.0761, 0.0000,  0.4976,  0.0000,  -0.4200,  0.0000},
                {     0.0216, 0.0000,  -0.1060, 0.0000,   0.0895,  0.0000},
                {     0.0254, 0.0000,  -0.1211, 0.0000,   0.1022,  0.0000},
                {    -0.2989, 0.0000,   1.3804, 0.0000,  -1.1650,  0.0000},
                {    -3.1873, 0.2010,  14.6890, 0.9266, -12.3974, -0.7820},
                {    -7.8468, 0.5320,  36.0910, 2.4469, -30.4606, -2.0652},
                {     0.0216, 0.0000,  -0.0988, 0.0000,   0.0834,  0.0000},
                {    -0.3384, 0.0000,   1.5433, 0.0000,  -1.3025,  0.0000},
                {     0.0179, 0.0000,  -0.0813, 0.0000,   0.0686,  0.0000},
                {    -0.0244, 0.0000,   0.1082, 0.0000,  -0.0913,  0.0000},
                {     0.0470, 0.0000,  -0.2004, 0.0000,   0.1692,  0.0000},
                {    -0.7341, 0.0000,   3.1240, 0.0000,  -2.6367,  0.0000},
                {    -0.0526, 0.0000,   0.2235, 0.0000,  -0.1886,  0.0000},
                {    -0.0508, 0.0000,   0.2073, 0.0000,  -0.1749,  0.0000},
                {     0.0498, 0.0000,  -0.1312, 0.0000,   0.1107,  0.0000},
                {     0.1006, 0.0000,  -0.2640, 0.0000,   0.2228,  0.0000},
                {     0.0395, 0.0000,  -0.0968, 0.0000,   0.0817,  0.0000},
                {     0.0470, 0.0000,  -0.1099, 0.0000,   0.0927,  0.0000},
                {     0.1767, 0.0000,  -0.4115, 0.0000,   0.3473,  0.0000},
                {     0.4352, 0.0000,  -1.0093, 0.0000,   0.8519,  0.0000},
                {     0.5339, 0.0000,  -1.2224, 0.0000,   1.0317,  0.0000},
                {    -8.4046, 0.2500,  19.1647, 0.5701, -16.1749, -0.4811},
                {     0.5443, 0.0000,  -1.2360, 0.0000,   1.0432,  0.0000},
                {     0.0470, 0.0000,  -0.1000, 0.0000,   0.0844,  0.0000},
                {    -0.0555, 0.0000,   0.1169, 0.0000,  -0.0987,  0.0000},
                {     0.1175, 0.0000,  -0.2332, 0.0000,   0.1968,  0.0000},
                {    -1.8236, 0.0000,   3.6018, 0.0000,  -3.0399,  0.0000},
                {     0.1316, 0.0000,  -0.2587, 0.0000,   0.2183,  0.0000},
                {     0.0179, 0.0000,  -0.0344, 0.0000,   0.0290,  0.0000},
                {    -0.0855, 0.0000,   0.1542, 0.0000,  -0.1302,  0.0000},
                {    -0.0573, 0.0000,   0.0395, 0.0000,  -0.0333,  0.0000},
                {     0.0329, 0.0000,  -0.0173, 0.0000,   0.0146,  0.0000},
                {    -1.8847, 0.0000,   0.9726, 0.0000,  -0.8209,  0.0000},
                {     0.2510, 0.0000,  -0.0910, 0.0000,   0.0768,  0.0000},
                {     1.1703, 0.0000,  -0.4135, 0.0000,   0.3490,  0.0000},
                {   -49.7174, 0.4330,  17.1056, 0.1490, -14.4370, -0.1257},
                {    -0.1936, 0.0000,   0.0666, 0.0000,  -0.0562,  0.0000},
                {     0.0489, 0.0000,  -0.0154, 0.0000,   0.0130,  0.0000},
                {    -0.5471, 0.0000,   0.1670, 0.0000,  -0.1409,  0.0000},
                {     0.0367, 0.0000,  -0.0108, 0.0000,   0.0092,  0.0000},
                {    -0.0451, 0.0000,   0.0082, 0.0000,  -0.0069,  0.0000},
                {     0.0921, 0.0000,  -0.0167, 0.0000,   0.0141,  0.0000},
                {     0.8281, 0.0000,  -0.1425, 0.0000,   0.1202,  0.0000},
                {   -15.8887, 0.1530,   2.7332, 0.0267,  -2.3068, -0.0222},
                {    -0.1382, 0.0000,   0.0225, 0.0000,  -0.0190,  0.0000},
                {     0.0348, 0.0000,  -0.0053, 0.0000,   0.0045,  0.0000},
                {    -0.1372, 0.0000,  -0.0079, 0.0000,   0.0066,  0.0000},
                {     0.4211, 0.0000,  -0.0203, 0.0000,   0.0171,  0.0000},
                {    -0.0404, 0.0000,   0.0008, 0.0000,  -0.0007,  0.0000},
                {     7.8998, 0.0000,   0.1460, 0.0000,  -0.1232,  0.0000},
                { -1617.2681, 0.0000, -14.9471, 0.0000,  12.6153,  0.0000}
            };
        }
        #endregion

        /// <summary>
        /// Constructor for calculation the zonal tides effect on Earth rotation.
        /// </summary>
        /// <param name="choise">Gives the component wihich has to be calculated.</param>
        public ZonalTides(Components Component)
        {
            this.Component = Component;

            switch (this.Component)
            {
                case Components.UT1: this.BaseUnit = BaseUnits.UT1; break;
                case Components.LOD: this.BaseUnit = BaseUnits.LOD; break;
                case Components.Omega: this.BaseUnit = BaseUnits.Omega; break;
                case Components.dCC: this.BaseUnit = BaseUnits.dCC; break;
            }
            _Unit = this.BaseUnit;
        }

        public double rg_zont2(double mjd)    //       SUBROUTINE RG_ZONT2 ( T, DUT, DLOD, DOMEGA )
        {
            /*  - - - - - - - - - - -
             *   R G _ Z O N T 2
             *  - - - - - - - - - - -

             *  This routine is part of the International Earth Rotation and
             *  Reference Systems Service (IERS) Conventions software collection.

             *  This subroutine evaluates the effects of zonal Earth tides on the
             *  rotation of the Earth.  The model used is a combination of Yoder
             *  et al. (1981) elastic body tide, Wahr and Bergen (1986) inelastic
             *  body tide, and Kantha et al. (1998) ocean tide models
             *  as recommended by the IERS Conventions (2010).  Refer to
             *  Chapter 8 pp. xx - xx.  The latest version of the model is located
             *  at http://tai.bipm.org/iers/convupdt/convupdt_c8.html.

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
             *     T           d      TT, Julian centuries since J2000 (Note 1)

             *  Returned:
             *     DUT         d      Effect on UT1 (Note 2)
             *     DLOD        d      Effect on excess length of day (LOD) (Note 3)
             *     DOMEGA      d      Effect on rotational speed (Note 4)

             *  Notes:

             *  1) Though T is strictly TDB, it is usually more convenient to use
             *     TT, which makes no significant difference.  Julian centuries since
             *     J2000 is (JD - 2451545.0)/36525.

             *  2) The expression used is as adopted in IERS Conventions (2010).
             *     DUT is expressed in seconds and is double precision.

             *  3) The expression used is as adopted in IERS Conventions (2010).
             *     DLOD is the excess in LOD and is expressed in seconds per day
             *     and is double precision.  The phrase 'per day' is generally
             *     understood, so it has been omitted commonly in speech and
             *     literature.
             *     See: Stephenson, F. R., Morrison, L. V., Whitrow, G. J., 1984,
             *     "Long-Term Changes in the Rotation of the Earth: 700 B. C. to
             *     A. D. 1980 [and Discussion]", Phil. Trans. Roy. Soc. of London.
             *     Series A, 313, pp. 47 - 70.

             *  4) The expression used is as adopted in IERS Conventions (2010).
             *     Rotational speed is expressed in radians per second and is
             *     double precision.

             *  Called:
             *     FUNDARG      Computation of the fundamental lunisolar arguments

             *  Test case:
             *     given input: T = .07995893223819302 Julian centuries since J2000
             *                  (MJD = 54465)
             *     expected output: DUT    =  7.983287678576557467E-002 seconds
             *                      DLOD   =  5.035303035410713729E-005 seconds / day
             *                      DOMEGA = -4.249711616463017E-014 radians / second

             *  References:

             *     Yoder, C. F., Williams, J. G., and Parke, M. E., (1981),
             *     "Tidal Variations of Earth Rotation," J. Geophys. Res., 86,
             *     pp. 881 - 891.

             *     Wahr, J. and Bergen, Z., (1986), "The effects of mantle
             *     anelasticity on nutations, Earth tides, and tidal variations
             *     in rotation rate," Geophys. J. Roy. astr. Soc., 87, pp. 633 - 668.

             *     Kantha, L. H., Stewart, J. S., and Desai, S. D., (1998), "Long-
             *     period lunar fortnightly and monthly ocean tides," J. Geophys.
             *     Res., 103, pp. 12639 - 12647.

             *     Gross, R. S., (2009), "Ocean tidal effects on Earth rotation,"
             *     J. Geodyn., 48(3-5), pp. 219 - 225.

             *     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
             *     IERS Technical Note No. xx, BKG (to be issued 2010)

             *  Revisions:
             *  2008 January 18 B.E. Stetzler  Initial changes to header
             *               and used 2PI instead of PI as parameter
             *  2008 January 25 B.E. Stetzler Additional changes to header
             *  2008 February 21 B.E. Stetzler Definition of (excess) LOD clarified
             *  2008 March   12 B.E. Stetzler Applied changes to wording of notes.
             *  2008 March   14 B.E. Stetzler Further changes applied to code.
             *  2008 April   03 B.E. Stetzler Provided example test case
             *  2009 February 11 B.E. Stetzler Updated test case due to changes made
             *                                 to FUNDARG.F subroutine
             *  2009 April   10 B.E. Stetzler DLOD corrected to say it is expressed
             *                                in seconds per day
             *  2009 May     04 B.E. Stetzler Code formatting changes based on
             *                                client recommendations
             *  2009 May     07 B.E. Stetzler Updated test case due to above changes
             *  2010 February 19 B.E. Stetzler Replaced Conventions 2003 recommended
             *                                 model with Conventions 2010 model
             *  2010 February 22 B.E. Stetzler Provided example test case
             *  2010 February 23 B.E. Stetzler Updated values to two decimal places
             *  2010 February 23 B.E. Stetzler Split fundamental arguments and
             *                                 coefficients for four decimal place
             *                                 precision
             *  2010 February 25 B.E. Stetzler Recalculation of fundamental arguments
             *  2010 March    01 B.E. Stetzler Updated table values to four decimal
             *                                 places and double precision
             *  2010 March    12 B.E. Stetzler Applied changes to wording of notes.
             *  2010 March    22 B.E. Stetzler Corrected DOMEGA output for test case
             * ----------------------------------------------------------------------- */

            /*   Computation of fundamental arguments */
            double L = 0.0, F = 0.0, D = 0.0, LP = 0.0, OM = 0.0;

            // compute julian century from J2000.0   
            mjd = (mjd - 51544.5d) / 36525.0d;

            CommonFunctions.fundarg(mjd, ref L, ref LP, ref F, ref D, ref OM);        //       CALL FUNDARG(T,L,LP,F,D,OM)

            /*  Set initial values to zero. */
            double dut = 0.0;       //       DUT    = 0.0D0
            double dlod = 0.0;      //       DLOD   = 0.0D0
            double domega = 0.0;        //       DOMEGA = 0.0D0

            /*  Sum zonal tide terms. */
            double ARG = 0.0;
            for (int i = 0; i < Table.NZONT; i++)     //       DO 10 I = 1, NZONT 
            {
                /*     Formation of multiples of arguments. */
                /*
                 ARG =      MOD ( DBLE ( NFUND(1,I) ) * L 
             .       +            DBLE ( NFUND(2,I) ) * LP 
             .       +            DBLE ( NFUND(3,I) ) * F 
             .       +            DBLE ( NFUND(4,I) ) * D 
             .       +            DBLE ( NFUND(5,I) ) * OM, D2PI ) */
                ARG = Table.nfund[i, 0] * L +
                      Table.nfund[i, 1] * LP +
                      Table.nfund[i, 2] * F +
                      Table.nfund[i, 3] * D +
                      Table.nfund[i, 4] * OM % Constants.TWOPI;

                //          IF (ARG.LT.0D0) ARG = ARG + D2PI

                if (ARG < 0.0)
                    ARG += Constants.TWOPI;

                /*     Evaluate zonal tidal terms. */

                dut = dut + Table.tide[i, 0] * Math.Sin(ARG) + Table.tide[i, 1] * Math.Cos(ARG);        //          DUT    = DUT    + TIDE(1,I) *DSIN(ARG) + TIDE(2,I) *DCOS(ARG)
                dlod = dlod + Table.tide[i, 2] * Math.Cos(ARG) + Table.tide[i, 3] * Math.Sin(ARG);      //          DLOD   = DLOD   + TIDE(3,I) *DCOS(ARG) + TIDE(4,I) *DSIN(ARG)
                domega = domega + Table.tide[i, 4] * Math.Cos(ARG) + Table.tide[i, 5] * Math.Sin(ARG);  //          DOMEGA = DOMEGA + TIDE(5,I) *DCOS(ARG) + TIDE(6,I) *DSIN(ARG)

            }       // 10    CONTINUE
            double result = double.NaN;

            /*  Rescale corrections so that they are in units involving seconds. */
            switch (Component)
            {
                case Components.UT1: // Correction to remove from UT1 in s
                    result = dut * 1.0e-4;
                    break;
                case Components.LOD: // Correction to remove from LOD in s
                    result = dlod * 1.0e-5;
                    break;
                case Components.Omega: // Correction to remove from Om in rad/s
                    result = domega * 1.0e-14;
                    break;
                case Components.dCC: //	Correction to remove from dC/C in rad
                    result = (domega * 1.0e-14) / Constants.Omega;
                    break;
            }

            return result * dConversion;
        }
    }
}
