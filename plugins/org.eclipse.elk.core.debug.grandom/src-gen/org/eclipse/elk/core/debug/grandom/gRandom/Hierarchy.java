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
 * A representation of the model object '<em><b>Hierarchy</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getLevels <em>Levels</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getEdges <em>Edges</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getNumHierarchNodes <em>Num Hierarch Nodes</em>}</li>
 *   <li>{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getCrossHierarchRel <em>Cross Hierarch Rel</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getHierarchy()
 * @model
 * @generated
 */
public interface Hierarchy extends EObject
{
  /**
   * Returns the value of the '<em><b>Levels</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Levels</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Levels</em>' containment reference.
   * @see #setLevels(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getHierarchy_Levels()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getLevels();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getLevels <em>Levels</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Levels</em>' containment reference.
   * @see #getLevels()
   * @generated
   */
  void setLevels(DoubleQuantity value);

  /**
   * Returns the value of the '<em><b>Edges</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Edges</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Edges</em>' containment reference.
   * @see #setEdges(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getHierarchy_Edges()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getEdges();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getEdges <em>Edges</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Edges</em>' containment reference.
   * @see #getEdges()
   * @generated
   */
  void setEdges(DoubleQuantity value);

  /**
   * Returns the value of the '<em><b>Num Hierarch Nodes</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Num Hierarch Nodes</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Num Hierarch Nodes</em>' containment reference.
   * @see #setNumHierarchNodes(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getHierarchy_NumHierarchNodes()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getNumHierarchNodes();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getNumHierarchNodes <em>Num Hierarch Nodes</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Num Hierarch Nodes</em>' containment reference.
   * @see #getNumHierarchNodes()
   * @generated
   */
  void setNumHierarchNodes(DoubleQuantity value);

  /**
   * Returns the value of the '<em><b>Cross Hierarch Rel</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Cross Hierarch Rel</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Cross Hierarch Rel</em>' containment reference.
   * @see #setCrossHierarchRel(DoubleQuantity)
   * @see org.eclipse.elk.core.debug.grandom.gRandom.GRandomPackage#getHierarchy_CrossHierarchRel()
   * @model containment="true"
   * @generated
   */
  DoubleQuantity getCrossHierarchRel();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.debug.grandom.gRandom.Hierarchy#getCrossHierarchRel <em>Cross Hierarch Rel</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Cross Hierarch Rel</em>' containment reference.
   * @see #getCrossHierarchRel()
   * @generated
   */
  void setCrossHierarchRel(DoubleQuantity value);

} // Hierarchy
