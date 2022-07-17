using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Elliptic_Curves
{
    /// <summary>
    /// An elliptic curve over Fp is the solution set to an equation of the form
    /// y^2 = x^3 + ax + b
    /// Where a,b in Fp satisfy 4a^3 + 27b^2 != 0
    /// This set also includes a point at infinity denoted O
    /// 
    /// We use E[a,b] to denote an elliptic curve y^2 = x^3 + ax + b
    /// 
    /// The points of the curve E[a,b] over Fp consists of all points (x,y) in Fp that satisfy the equation AND the point at infinity
    /// 
    /// If (x,y) is on E[a,b] then trivially (x, -y) is also
    /// 
    /// The condition 4a^3 + 27b^2 != 0 ensures that no repeated roots occur
    /// This enables us to define the points on the curve as a group
    /// In essence, we want to "add" two points, but this is done by taking tangents
    /// Tangents are not well defined when we have repeated roots
    /// I can highlight this clearer in person but see the the "Add" function
    /// </summary>
    internal class EllipticCurve
    {
        int discriminant, noOfSquares, a, b, p;
        List<Tuple<int, int>> pointsOnCurve = new List<Tuple<int, int>>();

        //Over R
        public EllipticCurve(int a, int b)
        {

        }
        //Over some field Fp, where p=mod, p prime
        public EllipticCurve(int a, int b, int mod)
        {
            this.a = a;
            this.b = b;
            p = mod;

            if (mod >= 5)
            {
                a = a % mod;
                b = b % mod;
                discriminant = (int)(4 * Math.Pow(a, 2) + 27 * Math.Pow(b, 3));
                if (discriminant != 0)
                {
                    noOfSquares = (int)((mod - 1.0) / 2) + 1;
                    generatePoints();
                }
                else
                {
                    throw new Exception();
                }
            }
        }
        private void generatePoints()
        {
            int x = 0;

            //We only need to find points where there are square roots
            //NOTE: point at infinity is assumed to exist and does not get added to list
            while (pointsOnCurve.Count<noOfSquares*2)
            {
                //Elliptic Curve formula
                int ySquared = (int)(Math.Pow(x, 3) + a * x + b) % p;

                //Use Eulers Criterion to see if square exists in Fp
                if (EulersCriterion(ySquared))
                {
                    int y = findRoot(ySquared);
                    pointsOnCurve.Add(new Tuple<int, int>(x, y));
                    //Add (x, -y)
                    //C#'s mod function does not handle negatives the way we want
                    //e.g. -1 % 19 == -1, but we want -1 % 19 == 18
                    //Since y<p we can get away with the following arithmetic instead
                    int yInv = -y + p;
                    pointsOnCurve.Add(new Tuple<int, int>(x, yInv));
                }
                x++;
            }
        }
        /// <summary>
        /// Famous theorem for finite fields
        /// A square exists for "a" iff a^[(p-1)/2] = 1 (mod p)
        /// </summary>
        /// <param name="a"></param>
        /// <returns>True if criterion holds, False otherwise</returns>
        public bool EulersCriterion(int a)
        {
            int power = (int)((p - 1.0) / 2);
            //powers can get big easily so ulong is needed
            ulong aToPower = (ulong)Math.Pow(a, power);
            if (aToPower%(ulong)p ==1)
                return true;
            return false;
        }
        public int findRoot(int ySquared)
        {
            for(int y = 0; y<p; y++)
            {
                //Find the y which is equal to ySquared
                if (Math.Pow(y, 2)%p == ySquared)
                    return y;
            }
            return 0;
        }
        static void Main(string[] args)
        {
            // Display the number of command line arguments.
            EllipticCurve curve = new EllipticCurve(1, 1, 19);
        }
    }
}
