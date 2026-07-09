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
 * A representation of the model object '<em><b>Edges</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isDensity <em>Density</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isTotal <em>Total</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isRelative <em>Relative</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isOutbound <em>Outbound</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#getNEdges <em>NEdges</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isLabels <em>Labels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isSelfLoops <em>Self Loops</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges()
 * @model
 * @generated
 */
public interface Edges extends EObject
{
  /**
   * Returns the value of the '<em><b>Density</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Density</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Density</em>' attribute.
   * @see #setDensity(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_Density()
   * @model
   * @generated
   */
  boolean isDensity();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isDensity <em>Density</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Density</em>' attribute.
   * @see #isDensity()
   * @generated
   */
  void setDensity(boolean value);

  /**
   * Returns the value of the '<em><b>Total</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Total</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Total</em>' attribute.
   * @see #setTotal(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_Total()
   * @model
   * @generated
   */
  boolean isTotal();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isTotal <em>Total</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Total</em>' attribute.
   * @see #isTotal()
   * @generated
   */
  void setTotal(boolean value);

  /**
   * Returns the value of the '<em><b>Relative</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Relative</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Relative</em>' attribute.
   * @see #setRelative(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_Relative()
   * @model
   * @generated
   */
  boolean isRelative();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isRelative <em>Relative</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Relative</em>' attribute.
   * @see #isRelative()
   * @generated
   */
  void setRelative(boolean value);

  /**
   * Returns the value of the '<em><b>Outbound</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Outbound</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Outbound</em>' attribute.
   * @see #setOutbound(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_Outbound()
   * @model
   * @generated
   */
  boolean isOutbound();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isOutbound <em>Outbound</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Outbound</em>' attribute.
   * @see #isOutbound()
   * @generated
   */
  void setOutbound(boolean value);

  /**
   * Returns the value of the '<em><b>NEdges</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>NEdges</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>NEdges</em>' containment reference.
   * @see #setNEdges(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_NEdges()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getNEdges();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#getNEdges <em>NEdges</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>NEdges</em>' containment reference.
   * @see #getNEdges()
   * @generated
   */
  void setNEdges(DoubleQuantity value);

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
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_Labels()
   * @model
   * @generated
   */
  boolean isLabels();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isLabels <em>Labels</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Labels</em>' attribute.
   * @see #isLabels()
   * @generated
   */
  void setLabels(boolean value);

  /**
   * Returns the value of the '<em><b>Self Loops</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Self Loops</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Self Loops</em>' attribute.
   * @see #setSelfLoops(boolean)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getEdges_SelfLoops()
   * @model
   * @generated
   */
  boolean isSelfLoops();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Edges#isSelfLoops <em>Self Loops</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Self Loops</em>' attribute.
   * @see #isSelfLoops()
   * @generated
   */
  void setSelfLoops(boolean value);

} // Edges
