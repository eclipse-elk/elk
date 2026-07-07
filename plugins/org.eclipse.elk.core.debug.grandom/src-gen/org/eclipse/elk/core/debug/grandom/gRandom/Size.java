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
 * A representation of the model object '<em><b>Size</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Size#getHeight <em>Height</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Size#getWidth <em>Width</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getSize()
 * @model
 * @generated
 */
public interface Size extends EObject
{
  /**
   * Returns the value of the '<em><b>Height</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Height</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Height</em>' containment reference.
   * @see #setHeight(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getSize_Height()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getHeight();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Size#getHeight <em>Height</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Height</em>' containment reference.
   * @see #getHeight()
   * @generated
   */
  void setHeight(DoubleQuantity value);

  /**
   * Returns the value of the '<em><b>Width</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Width</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Width</em>' containment reference.
   * @see #setWidth(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getSize_Width()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getWidth();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Size#getWidth <em>Width</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Width</em>' containment reference.
   * @see #getWidth()
   * @generated
   */
  void setWidth(DoubleQuantity value);

} // Size
