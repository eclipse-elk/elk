/**
 * ******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.core.debug.grandom.gRandom;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Double Quantity</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getQuant <em>Quant</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMin <em>Min</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isMinMax <em>Min Max</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMax <em>Max</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMean <em>Mean</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isGaussian <em>Gaussian</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getStddv <em>Stddv</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity()
 * @model
 * @generated
 */
public interface DoubleQuantity extends EObject
{
  /**
   * Returns the value of the '<em><b>Quant</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Quant</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Quant</em>' attribute.
   * @see #setQuant(Double)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_Quant()
   * @model
   * @generated
   */
  Double getQuant();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getQuant <em>Quant</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Quant</em>' attribute.
   * @see #getQuant()
   * @generated
   */
  void setQuant(Double value);

  /**
   * Returns the value of the '<em><b>Min</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Min</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Min</em>' attribute.
   * @see #setMin(Double)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_Min()
   * @model
   * @generated
   */
  Double getMin();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMin <em>Min</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Min</em>' attribute.
   * @see #getMin()
   * @generated
   */
  void setMin(Double value);

  /**
   * Returns the value of the '<em><b>Min Max</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Min Max</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Min Max</em>' attribute.
   * @see #setMinMax(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_MinMax()
   * @model
   * @generated
   */
  boolean isMinMax();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isMinMax <em>Min Max</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Min Max</em>' attribute.
   * @see #isMinMax()
   * @generated
   */
  void setMinMax(boolean value);

  /**
   * Returns the value of the '<em><b>Max</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Max</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Max</em>' attribute.
   * @see #setMax(Double)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_Max()
   * @model
   * @generated
   */
  Double getMax();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMax <em>Max</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Max</em>' attribute.
   * @see #getMax()
   * @generated
   */
  void setMax(Double value);

  /**
   * Returns the value of the '<em><b>Mean</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Mean</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Mean</em>' attribute.
   * @see #setMean(Double)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_Mean()
   * @model
   * @generated
   */
  Double getMean();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getMean <em>Mean</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Mean</em>' attribute.
   * @see #getMean()
   * @generated
   */
  void setMean(Double value);

  /**
   * Returns the value of the '<em><b>Gaussian</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Gaussian</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Gaussian</em>' attribute.
   * @see #setGaussian(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_Gaussian()
   * @model
   * @generated
   */
  boolean isGaussian();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#isGaussian <em>Gaussian</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Gaussian</em>' attribute.
   * @see #isGaussian()
   * @generated
   */
  void setGaussian(boolean value);

  /**
   * Returns the value of the '<em><b>Stddv</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Stddv</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Stddv</em>' attribute.
   * @see #setStddv(Double)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getDoubleQuantity_Stddv()
   * @model
   * @generated
   */
  Double getStddv();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.DoubleQuantity#getStddv <em>Stddv</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Stddv</em>' attribute.
   * @see #getStddv()
   * @generated
   */
  void setStddv(Double value);

} // DoubleQuantity
