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

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Rand Graph</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.RandGraph#getConfigs <em>Configs</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getRandGraph()
 * @model
 * @generated
 */
public interface RandGraph extends EObject
{
  /**
   * Returns the value of the '<em><b>Configs</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.core.debug.grandom.gRandom.Configuration}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Configs</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Configs</em>' containment reference list.
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getRandGraph_Configs()
   * @model containment="true"
   * @generated
   */
  EList<Configuration> getConfigs();

} // RandGraph
