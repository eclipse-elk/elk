/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom.impl;

import org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity;
import org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage;

import org.eclipse.emf.common.notify.Notification;

import org.eclipse.emf.ecore.EClass;

import org.eclipse.emf.ecore.impl.ENotificationImpl;
import org.eclipse.emf.ecore.impl.MinimalEObjectImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model object '<em><b>Double Quantity</b></em>'.
 * <!-- end-user-doc -->
 * <p>
 * The following features are implemented:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#getQuant <em>Quant</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#getMin <em>Min</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#isMinMax <em>Min Max</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#getMax <em>Max</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#getMean <em>Mean</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#isGaussian <em>Gaussian</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.impl.DoubleQuantityImpl#getStddv <em>Stddv</em>}</li>
 * </ul>
 *
 * @generated
 */
public class DoubleQuantityImpl extends MinimalEObjectImpl.Container implements DoubleQuantity
{
  /**
   * The default value of the '{@link #getQuant() <em>Quant</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getQuant()
   * @generated
   * @ordered
   */
  protected static final Double QUANT_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getQuant() <em>Quant</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getQuant()
   * @generated
   * @ordered
   */
  protected Double quant = QUANT_EDEFAULT;

  /**
   * The default value of the '{@link #getMin() <em>Min</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMin()
   * @generated
   * @ordered
   */
  protected static final Double MIN_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getMin() <em>Min</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMin()
   * @generated
   * @ordered
   */
  protected Double min = MIN_EDEFAULT;

  /**
   * The default value of the '{@link #isMinMax() <em>Min Max</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isMinMax()
   * @generated
   * @ordered
   */
  protected static final boolean MIN_MAX_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isMinMax() <em>Min Max</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isMinMax()
   * @generated
   * @ordered
   */
  protected boolean minMax = MIN_MAX_EDEFAULT;

  /**
   * The default value of the '{@link #getMax() <em>Max</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMax()
   * @generated
   * @ordered
   */
  protected static final Double MAX_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getMax() <em>Max</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMax()
   * @generated
   * @ordered
   */
  protected Double max = MAX_EDEFAULT;

  /**
   * The default value of the '{@link #getMean() <em>Mean</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMean()
   * @generated
   * @ordered
   */
  protected static final Double MEAN_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getMean() <em>Mean</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getMean()
   * @generated
   * @ordered
   */
  protected Double mean = MEAN_EDEFAULT;

  /**
   * The default value of the '{@link #isGaussian() <em>Gaussian</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isGaussian()
   * @generated
   * @ordered
   */
  protected static final boolean GAUSSIAN_EDEFAULT = false;

  /**
   * The cached value of the '{@link #isGaussian() <em>Gaussian</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #isGaussian()
   * @generated
   * @ordered
   */
  protected boolean gaussian = GAUSSIAN_EDEFAULT;

  /**
   * The default value of the '{@link #getStddv() <em>Stddv</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getStddv()
   * @generated
   * @ordered
   */
  protected static final Double STDDV_EDEFAULT = null;

  /**
   * The cached value of the '{@link #getStddv() <em>Stddv</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @see #getStddv()
   * @generated
   * @ordered
   */
  protected Double stddv = STDDV_EDEFAULT;

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected DoubleQuantityImpl()
  {
    super();
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  protected EClass eStaticClass()
  {
    return GRandomPackage.Literals.DOUBLE_QUANTITY;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Double getQuant()
  {
    return quant;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setQuant(Double newQuant)
  {
    Double oldQuant = quant;
    quant = newQuant;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__QUANT, oldQuant, quant));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Double getMin()
  {
    return min;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setMin(Double newMin)
  {
    Double oldMin = min;
    min = newMin;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__MIN, oldMin, min));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isMinMax()
  {
    return minMax;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setMinMax(boolean newMinMax)
  {
    boolean oldMinMax = minMax;
    minMax = newMinMax;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__MIN_MAX, oldMinMax, minMax));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Double getMax()
  {
    return max;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setMax(Double newMax)
  {
    Double oldMax = max;
    max = newMax;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__MAX, oldMax, max));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Double getMean()
  {
    return mean;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setMean(Double newMean)
  {
    Double oldMean = mean;
    mean = newMean;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__MEAN, oldMean, mean));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public boolean isGaussian()
  {
    return gaussian;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setGaussian(boolean newGaussian)
  {
    boolean oldGaussian = gaussian;
    gaussian = newGaussian;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__GAUSSIAN, oldGaussian, gaussian));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public Double getStddv()
  {
    return stddv;
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public void setStddv(Double newStddv)
  {
    Double oldStddv = stddv;
    stddv = newStddv;
    if (eNotificationRequired())
      eNotify(new ENotificationImpl(this, Notification.SET, GRandomPackage.DOUBLE_QUANTITY__STDDV, oldStddv, stddv));
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public Object eGet(int featureID, boolean resolve, boolean coreType)
  {
    switch (featureID)
    {
      case GRandomPackage.DOUBLE_QUANTITY__QUANT:
        return getQuant();
      case GRandomPackage.DOUBLE_QUANTITY__MIN:
        return getMin();
      case GRandomPackage.DOUBLE_QUANTITY__MIN_MAX:
        return isMinMax();
      case GRandomPackage.DOUBLE_QUANTITY__MAX:
        return getMax();
      case GRandomPackage.DOUBLE_QUANTITY__MEAN:
        return getMean();
      case GRandomPackage.DOUBLE_QUANTITY__GAUSSIAN:
        return isGaussian();
      case GRandomPackage.DOUBLE_QUANTITY__STDDV:
        return getStddv();
    }
    return super.eGet(featureID, resolve, coreType);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public void eSet(int featureID, Object newValue)
  {
    switch (featureID)
    {
      case GRandomPackage.DOUBLE_QUANTITY__QUANT:
        setQuant((Double)newValue);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MIN:
        setMin((Double)newValue);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MIN_MAX:
        setMinMax((Boolean)newValue);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MAX:
        setMax((Double)newValue);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MEAN:
        setMean((Double)newValue);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__GAUSSIAN:
        setGaussian((Boolean)newValue);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__STDDV:
        setStddv((Double)newValue);
        return;
    }
    super.eSet(featureID, newValue);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public void eUnset(int featureID)
  {
    switch (featureID)
    {
      case GRandomPackage.DOUBLE_QUANTITY__QUANT:
        setQuant(QUANT_EDEFAULT);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MIN:
        setMin(MIN_EDEFAULT);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MIN_MAX:
        setMinMax(MIN_MAX_EDEFAULT);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MAX:
        setMax(MAX_EDEFAULT);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__MEAN:
        setMean(MEAN_EDEFAULT);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__GAUSSIAN:
        setGaussian(GAUSSIAN_EDEFAULT);
        return;
      case GRandomPackage.DOUBLE_QUANTITY__STDDV:
        setStddv(STDDV_EDEFAULT);
        return;
    }
    super.eUnset(featureID);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public boolean eIsSet(int featureID)
  {
    switch (featureID)
    {
      case GRandomPackage.DOUBLE_QUANTITY__QUANT:
        return QUANT_EDEFAULT == null ? quant != null : !QUANT_EDEFAULT.equals(quant);
      case GRandomPackage.DOUBLE_QUANTITY__MIN:
        return MIN_EDEFAULT == null ? min != null : !MIN_EDEFAULT.equals(min);
      case GRandomPackage.DOUBLE_QUANTITY__MIN_MAX:
        return minMax != MIN_MAX_EDEFAULT;
      case GRandomPackage.DOUBLE_QUANTITY__MAX:
        return MAX_EDEFAULT == null ? max != null : !MAX_EDEFAULT.equals(max);
      case GRandomPackage.DOUBLE_QUANTITY__MEAN:
        return MEAN_EDEFAULT == null ? mean != null : !MEAN_EDEFAULT.equals(mean);
      case GRandomPackage.DOUBLE_QUANTITY__GAUSSIAN:
        return gaussian != GAUSSIAN_EDEFAULT;
      case GRandomPackage.DOUBLE_QUANTITY__STDDV:
        return STDDV_EDEFAULT == null ? stddv != null : !STDDV_EDEFAULT.equals(stddv);
    }
    return super.eIsSet(featureID);
  }

  /**
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  @Override
  public String toString()
  {
    if (eIsProxy()) return super.toString();

    StringBuilder result = new StringBuilder(super.toString());
    result.append(" (quant: ");
    result.append(quant);
    result.append(", min: ");
    result.append(min);
    result.append(", minMax: ");
    result.append(minMax);
    result.append(", max: ");
    result.append(max);
    result.append(", mean: ");
    result.append(mean);
    result.append(", gaussian: ");
    result.append(gaussian);
    result.append(", stddv: ");
    result.append(stddv);
    result.append(')');
    return result.toString();
  }

} //DoubleQuantityImpl
