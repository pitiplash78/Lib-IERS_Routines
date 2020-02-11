using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IERS_Routines
{
    public static class  Constants
    {
        /// <summary>
        /// Constant number PI
        /// </summary>
        public static readonly double PI = 3.141592653589793238462643;
        
        /// <summary>
        /// Constant number 2xPI
        /// </summary>
        public static readonly double TWOPI = 6.283185307179586476925287;
                                              
        /// <summary>
        /// modified Julian date of J2000
        /// </summary>
        public static readonly double RMJD0 = 51544.5;

        /// <summary>
        /// Radians to seconds
        /// </summary>
        public static readonly double RAD2SEC = 13750.987083139758;

        /// <summary>
        /// Arcseconds to radians
        /// </summary>
        public static readonly double DAS2R = 4.848136811095359935899141e-6;

        /// <summary>
        /// Arcseconds in a full circle
        /// </summary>
        public static readonly double TURNAS = 1.296e6;// 1296000.0

        /// <summary>
        /// Mean angular velocity od the Earth in rad/s. 
        /// </summary>
        public static readonly double Omega = 7.292115e-5;
    }
}
