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
package org.eclipse.elk.core.debug.grandom.gRandom;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Flow</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getFlowType <em>Flow Type</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getSide <em>Side</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getAmount <em>Amount</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getFlow()
 * @model
 * @generated
 */
public interface Flow extends EObject
{
  /**
   * Returns the value of the '<em><b>Flow Type</b></em>' attribute.
   * The literals are from the enumeration {@link org.eclipse.elk.core.debug.grandom.gRandom.FlowType}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Flow Type</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Flow Type</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.FlowType
   * @see #setFlowType(FlowType)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getFlow_FlowType()
   * @model
   * @generated
   */
  FlowType getFlowType();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getFlowType <em>Flow Type</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Flow Type</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.FlowType
   * @see #getFlowType()
   * @generated
   */
  void setFlowType(FlowType value);

  /**
   * Returns the value of the '<em><b>Side</b></em>' attribute.
   * The literals are from the enumeration {@link org.eclipse.elk.core.debug.grandom.gRandom.Side}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Side</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Side</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Side
   * @see #setSide(Side)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getFlow_Side()
   * @model
   * @generated
   */
  Side getSide();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getSide <em>Side</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Side</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.Side
   * @see #getSide()
   * @generated
   */
  void setSide(Side value);

  /**
   * Returns the value of the '<em><b>Amount</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Amount</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Amount</em>' containment reference.
   * @see #setAmount(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getFlow_Amount()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getAmount();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Flow#getAmount <em>Amount</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Amount</em>' containment reference.
   * @see #getAmount()
   * @generated
   */
  void setAmount(DoubleQuantity value);

} // Flow
