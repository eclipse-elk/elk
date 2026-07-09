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

import org.eclipse.emf.ecore.EFactory;

/**
 * <!-- begin-user-doc -->
 * The <b>Factory</b> for the model.
 * It provides a create method for each non-abstract class of the model.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage
 * @generated
 */
public interface GRandomFactory extends EFactory
{
  /**
   * The singleton instance of the factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  GRandomFactory eINSTANCE = org.eclipse.elk.core.debug.grandom.gRandom.impl.GRandomFactoryImpl.init();

  /**
   * Returns a new object of class '<em>Rand Graph</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Rand Graph</em>'.
   * @generated
   */
  RandGraph createRandGraph();

  /**
   * Returns a new object of class '<em>Configuration</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Configuration</em>'.
   * @generated
   */
  Configuration createConfiguration();

  /**
   * Returns a new object of class '<em>Hierarchy</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Hierarchy</em>'.
   * @generated
   */
  Hierarchy createHierarchy();

  /**
   * Returns a new object of class '<em>Edges</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Edges</em>'.
   * @generated
   */
  Edges createEdges();

  /**
   * Returns a new object of class '<em>Nodes</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Nodes</em>'.
   * @generated
   */
  Nodes createNodes();

  /**
   * Returns a new object of class '<em>Size</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Size</em>'.
   * @generated
   */
  Size createSize();

  /**
   * Returns a new object of class '<em>Ports</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Ports</em>'.
   * @generated
   */
  Ports createPorts();

  /**
   * Returns a new object of class '<em>Flow</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Flow</em>'.
   * @generated
   */
  Flow createFlow();

  /**
   * Returns a new object of class '<em>Double Quantity</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Double Quantity</em>'.
   * @generated
   */
  DoubleQuantity createDoubleQuantity();

  /**
   * Returns the package supported by this factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the package supported by this factory.
   * @generated
   */
  GRandomPackage getGRandomPackage();

} //GRandomFactory
