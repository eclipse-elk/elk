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

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Ports</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#isLabels <em>Labels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getReUse <em>Re Use</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getSize <em>Size</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getConstraint <em>Constraint</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getFlow <em>Flow</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getPorts()
 * @model
 * @generated
 */
public interface Ports extends EObject
{
  /**
   * Returns the value of the '<em><b>Labels</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Labels</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Labels</em>' attribute.
   * @see #setLabels(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getPorts_Labels()
   * @model
   * @generated
   */
  boolean isLabels();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#isLabels <em>Labels</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Labels</em>' attribute.
   * @see #isLabels()
   * @generated
   */
  void setLabels(boolean value);

  /**
   * Returns the value of the '<em><b>Re Use</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Re Use</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Re Use</em>' containment reference.
   * @see #setReUse(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getPorts_ReUse()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getReUse();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getReUse <em>Re Use</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Re Use</em>' containment reference.
   * @see #getReUse()
   * @generated
   */
  void setReUse(DoubleQuantity value);

  /**
   * Returns the value of the '<em><b>Size</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Size</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Size</em>' containment reference.
   * @see #setSize(Size)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getPorts_Size()
   * @model containment="true"
   * @generated
   */
  Size getSize();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getSize <em>Size</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Size</em>' containment reference.
   * @see #getSize()
   * @generated
   */
  void setSize(Size value);

  /**
   * Returns the value of the '<em><b>Constraint</b></em>' attribute.
   * The literals are from the enumeration {@link org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Constraint</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Constraint</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType
   * @see #setConstraint(ConstraintType)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getPorts_Constraint()
   * @model
   * @generated
   */
  ConstraintType getConstraint();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Ports#getConstraint <em>Constraint</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Constraint</em>' attribute.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.ConstraintType
   * @see #getConstraint()
   * @generated
   */
  void setConstraint(ConstraintType value);

  /**
   * Returns the value of the '<em><b>Flow</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.core.debug.grandom.gRandom.Flow}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Flow</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Flow</em>' containment reference list.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getPorts_Flow()
   * @model containment="true"
   * @generated
   */
  EList<Flow> getFlow();

} // Ports
