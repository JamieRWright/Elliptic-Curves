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
    /// 
    /// This set also includes a point at infinity denoted O
    /// This point at infinity is to be thought of as being located far away along the y-direction.
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
    /// 
    /// 
    /// </summary>
    internal class EllipticCurve 
    {
        int discriminant, noOfSquares, ord, a, b, p;
        int[] inverses;
        List<Tuple<int, int>> pointsOnCurve = new List<Tuple<int, int>>();
        Tuple<int, int> generator = default;


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
                    ord = noOfSquares*2 + 1;
                    generatePoints();
                    generateInverses();

                    generator = getGenerator();
                }
                else
                {
                    throw new Exception();
                }
            }
        }


        public EllipticCurve(int a, int b, int mod, Tuple<int, int> inGenerator)
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
                    ord = noOfSquares * 2 + 1;
                    generatePoints();
                    generateInverses();

                    generator = inGenerator;
                }
                else
                {
                    throw new Exception();
                }
            }
        }


        private Tuple<int, int>? getGenerator()
        {
            foreach (Tuple<int, int> point in pointsOnCurve)
            {
                if (isGenerator(point))
                    return point;
            }
            //Generators don't always exist
            return null;
        }
        private void generateInverses()
        {
            bool isInverse = false;
            inverses = new int[p];
            inverses[0] = 0;

            for(int i = 1; i < p; i++)
            {
                int inv = 1;
                
                do
                {
                    isInverse = (i * inv) % p == 1;
                    if (isInverse)
                        inverses[i] = inv;
                    else
                        inv++;
                }
                while (inv < p && !isInverse);
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

        /// <summary>
        /// Division works differently in the world of modulo arithmetic so a special function is required
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public int div(int numerator, int denominator)
        {
            //First make any negative numbers positive via modding
            if (numerator < 0)
                numerator = numerator + p;
            if (denominator < 0)
                denominator = denominator + p;

            //Find the inverse for the given denominator and multiply with the numerator mod p
            int denominator_inv = inverses[denominator]; 

            //now find inverse of b
            return (numerator * denominator_inv) % p;
        }
        /// <summary>
        /// Generate binary operation on an elliptic curve over Fp which satisfy properties of an abelian group.
        /// Consider geometrically the addition of two points on an elliptic curve to be the following:
        /// First take a line through P1 and P2, this line will intersect the curve at a third pont.
        /// The reflection of this point around the x-axis (invert y-coord) also lies on the curve and is the point P1 + P2.
        /// </summary>
        /// <param name="P1"></param>
        /// <param name="P2"></param>
        /// <returns></returns>
        public Tuple<int, int> Add(Tuple<int, int> P1, Tuple<int, int> P2)
        {
            //This method makes me hate C#'s division of integers :(

            //Assign x1, y1, x2, y2 for simplicity
            int x1, y1, x2, y2, x3 = -1, y3 = -1, m, c;
            x1 = P1.Item1;
            y1 = P1.Item2;
            x2 = P2.Item1;
            y2 = P2.Item2;

            //Booleans to check if P1 or P2 are points at infinity
            bool P1_PAI = (x1 == p && y1 == p);
            bool P2_PAI = (x2 == p && y2 == p);

            //Case 1: x1 != x2
            if (x1 != x2)
            {
                //It is clear that the points P1, P2 satisfy both y = mx + c & E[a,b]
                //This is basic manipulation of the equation of a line between two distinct points
                m = div((y1 - y2),  (x1 - x2));
                c = div((x2 * y1 - x1 * y2), (x2 - x1));

                //Subsituting the sim equations we obtain that x^3 + m^2(x^2) ...
                //This is all we need as we know for cubic equations the roots add up to this m^2
                x3 = (int)(Math.Pow(m, 2) - x1 - x2) % p;

                //Sub into y=mx+c but flip the y
                y3 = -1 * (m * x3 + c) % p;
            }

            //Case 2: P1 = P2 & (y != 0) i.e. 2P1
            else if (x1 == x2 && y1 == y2 && y1 != 0)
            {
                //For this case we need to find a tangent at the point, so implicit differentation is needed
                //Differentiating E[a,b] and substituting in the points P1
                m = div((3 * (int)Math.Pow(x1,2) + a), (2 * y1));
                c = (y1 - m * x1) % p;
                if (c < 0)
                    c = c + p;

                x3 = (int)(Math.Pow(m, 2) - 2 * x1) % p;

                //Sub into y=mx+c but flip the y
                y3 = -1 * (m * x3 + c) % p;

            }

            //Case 3: P1 = (x1, y1), P2 = (x1, -y1) where y can equal 0
            else if (x1 == x2 && y1 == (-1*y2 + p))
            {
                x3 = p;
                y3 = p;
            }

            //Case 4: Where P1 or P2 is the point at infinity
            else if (P1_PAI || P2_PAI)
            {
                if (P1_PAI)
                    return P2;
                else
                    return P1;
            }

            if (x3 < 0)
                x3 = x3 + p;
            if (y3 < 0)
                y3 = y3 + p;
            return new Tuple<int, int>(x3, y3);
        }
        // Generate binary operation on an elliptic curve over Fp which satisfy properties of an abelian group
        // Explicitly 5 properties will hold:
        // 
        // 1: For P1, P2 in E, then P1 + P2 is in E (Closure under addition)
        // 2: The element O in E exists and holds such that O + P is in E for all P in E (Identity)
        // 3: P1 + (P2 + P3) = (P1 + P2) + P3 for all P1, P2, P3 in E (Assosiative Law)
        // 4: For all (x, y) in E, there exists -(x, y) in E such that (x, y) + -(x, y) = O (when (x,y) != O, this is simply (x, -y)) (Inverse)
        // 
        // We now have a group, but it is an abelian group as the following holds
        // 5: P1 + P2 = P2 + P1 for all P1, P2 in E 

        /// <summary>
        /// Adds the point P together k times
        /// </summary>
        /// <param name="k"></param>
        /// <param name="P"></param>
        /// <returns>Returns the point kP</returns>
        public Tuple<int, int> multiplyPoint(int k, Tuple<int, int> P)
        {
            bool isInverse = false;
            Tuple<int, int> kP = P;

            //if negative, invert at end
            if (k < 0)
            {
                k = Math.Abs(k);
                isInverse = true;
            }

            if (k == 0)
                return new Tuple<int, int>(p, p);

            else if (k == 1)
            {
                if (isInverse)
                    return new Tuple<int, int>(P.Item1, -P.Item2 + p);
                else 
                    return P;
            }
            else
            {
                for (int i = 2; i <= k; i++)
                {
                    kP = Add(kP, P);
                }

                if (isInverse)
                    return new Tuple<int, int>(kP.Item1, -kP.Item2 + p);
                else
                    return kP;
            }
        }

        /// <summary>
        /// Checks if a point is a generator
        /// </summary>
        /// <param name="P"></param>
        /// <returns>True if P is a generator, False otherwise</returns>
        public bool isGenerator(Tuple<int, int> P)
        {
            if (ordOfP(P) == ord)
                return true;
            return false;
        }

        /// <summary>
        /// Finds all the factors of a number
        /// </summary>
        /// <param name="number"></param>
        /// <returns>List containing the factors of number</returns>
        public List<int> Factor(int number)
        {
            var factors = new List<int>();
            int max = (int)Math.Sqrt(number);  // Round down

            for (int factor = 1; factor <= max; ++factor) // Test from 1 to the square root, or the int below it, inclusive.
            {
                if (number % factor == 0)
                {
                    factors.Add(factor);
                    if (factor != number / factor) // Don't add the square root twice!  Thanks Jon
                        factors.Add(number / factor);
                }
            }
            return factors;
        }


        /// <summary>
        /// Finds the order of a given point P
        /// </summary>
        /// <param name="P"></param>
        public int ordOfP(Tuple<int, int> P)
        {
            Tuple<int, int> kP = P;

            //Order of a point can only be a divisor of p
            foreach(int k in Factor(ord))
            {
                kP = multiplyPoint(k, P);
                if (kP.Item1 == p && kP.Item2 == p)
                    return k;
            }
            return -1;
        }


        /// <summary>
        /// Using a secret key from a sender and a published public key from a reciepient
        /// A message is encrypted using El Gamal
        /// </summary>
        /// <param name="secretKey"></param>
        /// <param name="publicKey"></param>
        /// <param name="M"></param>
        /// <returns>The pair of points (kP, M + kPI) is returned </returns>
        public Tuple<Tuple<int, int>, Tuple<int, int>> ElGamalSend(int secretKey, Tuple<int, int> publicKey, Tuple<int, int> M)
        {

            Tuple<int, int> kP = multiplyPoint(secretKey, generator);
            Tuple<int, int> M_Add_kPI = Add(M, multiplyPoint(secretKey, publicKey));

            return new Tuple<Tuple<int, int>, Tuple<int, int>>(kP, M_Add_kPI);
        }

        /// <summary>
        /// Decrypts the message sent by El Gamal using the pair of points recieved
        /// </summary>
        /// <param name="secretKey"></param>
        /// <param name="kP"></param>
        /// <param name="M_Add_kPI"></param>
        /// <returns>Message M</returns>
        public Tuple<int, int> ElGamalRecieve(int secretKey, Tuple<int, int> kP, Tuple<int, int> M_Add_kPI)
        {
            Tuple<int, int> M;
            Tuple<int, int> kaP = multiplyPoint(secretKey, kP);
            //Invert this point
            kaP = new Tuple<int, int>(kaP.Item1, -kaP.Item2);


            M = Add(M_Add_kPI, kaP); 
            return M;
        }


        static void Main(string[] args)
        {
            // Display the number of command line arguments.
            EllipticCurve curve = new EllipticCurve(1, 1, 19);
            Tuple<int, int> P1 = curve.pointsOnCurve[0];
            Tuple<int, int> P2 = curve.pointsOnCurve[8];
            Tuple<int, int> P1_Add_P2 = curve.Add(P1, P2);


            //EL GAMAL

            Tuple<int, int> publicKey = curve.pointsOnCurve[8];
            Tuple<int, int> M = new Tuple<int, int>(2, 7);

            Tuple<Tuple<int, int>, Tuple<int, int>> returnPair = curve.ElGamalSend(2, publicKey, M);
            Tuple<int, int> M_Recieved = curve.ElGamalRecieve(17, returnPair.Item1, returnPair.Item2);
            Console.WriteLine(P1_Add_P2);

        }
    }
}
