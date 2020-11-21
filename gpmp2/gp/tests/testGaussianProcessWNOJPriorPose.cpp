/**
*  @file testGaussianProcessWNOJPriorPose2.cpp
*  @author Jing Dong
**/

#include <CppUnitLite/TestHarness.h>

#include <gtsam/base/Matrix.h>
#include <gtsam/base/Testable.h>
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>

#include <gpmp2/gp/GaussianProcessWNOJPriorPose.h>

#include <iostream>

using namespace std;
using namespace gtsam;
using namespace gpmp2;


/* ************************************************************************** */
TEST(GaussianProcessWNOJPriorPose2, Factor) {

  const double delta_t = 0.1;
  Matrix Qc = 0.01 * Matrix::Identity(3,3);
  noiseModel::Gaussian::shared_ptr Qc_model = noiseModel::Gaussian::Covariance(Qc);
  Key key_pose1 = Symbol('x', 1), key_pose2 = Symbol('x', 2);
  Key key_vel1 = Symbol('v', 1), key_vel2 = Symbol('v', 2);
  Key key_acc1 = Symbol('a', 1), key_acc2 = Symbol('a', 2);
  GaussianProcessWNOJPriorPose2 factor(key_pose1, key_vel1, key_acc1, key_pose2, key_vel2, key_acc1, delta_t, Qc_model);
  Pose2 p1, p2;
  Vector3 v1, v2, a1, a2;
  Matrix actualH1, actualH2, actualH3, actualH4, actualH5, actualH6;
  Matrix expectH1, expectH2, expectH3, expectH4, expectH5, expectH6;
  Vector actual, expect;


  // test at origin
  p1 = Pose2(0, 0, 0);
  p2 = Pose2(0, 0, 0);
  v1 = (Vector3() << 0, 0, 0).finished();
  v2 = (Vector3() << 0, 0, 0).finished();
  a1 = (Vector3() << 0, 0, 0).finished();
  a2 = (Vector3() << 0, 0, 0).finished();
  actual = factor.evaluateError(p1, v1, a1, p2, v2, a2, actualH1, actualH2, actualH3, actualH4 ,actualH5, actualH6);
  expect = (Vector(6) << 0, 0, 0, 0, 0, 0).finished();

  // helper function
  boost::function<Pose2(const Pose2 &, const Vector3 &, const Vector3 &,
                        const Pose2 &, const Vector3 &, const Vector3 &)>
      f1 = [factor](const Pose2 &p1, const Vector3 &v1, const Vector3 &a1,
                  const Pose2 &p2, const Vector3 &v2, const Vector3 &a2) {
        return base.interpolatePose(p1, v1, a1, p2, v2, a2);
      };

  // Calculate derivates using numerical approximation
  expectH1 = numericalDerivative61(f1, p1, v1, a1, p2, v2, a2, 1e-6);
  expectH2 = numericalDerivative62(f1, p1, v1, a1, p2, v2, a2, 1e-6);
  expectH3 = numericalDerivative63(f1, p1, v1, a1, p2, v2, a2, 1e-6);
  expectH4 = numericalDerivative64(f1, p1, v1, a1, p2, v2, a2, 1e-6);
  expectH5 = numericalDerivative65(f1, p1, v1, a1, p2, v2, a2, 1e-6);
  expectH6 = numericalDerivative66(f1, p1, v1, a1, p2, v2, a2, 1e-6);

  EXPECT(assert_equal(expect, actual, 1e-6));
  EXPECT(assert_equal(expectH1, actualH1, 1e-8));
  EXPECT(assert_equal(expectH2, actualH2, 1e-8));
  EXPECT(assert_equal(expectH3, actualH3, 1e-8));
  EXPECT(assert_equal(expectH4, actualH4, 1e-8));
  EXPECT(assert_equal(expectH5, actualH5, 1e-8));
  EXPECT(assert_equal(expectH6, actualH6, 1e-8));


}

/* ************************************************************************** */
TEST(GaussianProcessWNOJPriorPose2, Optimization) {
  /**
   * A simple graph:
   *
   * p1   p2
   * |    |
   * x1   x2
   *  \  /
   * a1-gp-a2
   *  /  \
   * v1  v2
   *
   * p1 and p2 are pose prior factor to fix the poses, gp is the GP factor
   * that get correct velocity of v2
   */

  
}

/* ************************************************************************** */
/* main function */
int main() {
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
