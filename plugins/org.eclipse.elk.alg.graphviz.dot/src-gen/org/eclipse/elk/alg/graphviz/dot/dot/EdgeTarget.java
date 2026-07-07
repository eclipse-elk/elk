/**
 * ******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *  ******************************************************************************
 */
package org.eclipse.elk.alg.graphviz.dot.dot;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Edge Target</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget#getOperator <em>Operator</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget#getTargetSubgraph <em>Target Subgraph</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget#getTargetnode <em>Targetnode</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeTarget()
 * @model
 * @generated
 */
public interface EdgeTarget extends EObject
{
  /**
   * Returns the value of the '<em><b>Operator</b></em>' attribute.
   * The literals are from the enumeration {@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeOperator}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Operator</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Operator</em>' attribute.
   * @see org.eclipse.elk.alg.graphviz.dot.dot.EdgeOperator
   * @see #setOperator(EdgeOperator)
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeTarget_Operator()
   * @model
   * @generated
   */
  EdgeOperator getOperator();

  /**
   * Sets the value of the '{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget#getOperator <em>Operator</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Operator</em>' attribute.
   * @see org.eclipse.elk.alg.graphviz.dot.dot.EdgeOperator
   * @see #getOperator()
   * @generated
   */
  void setOperator(EdgeOperator value);

  /**
   * Returns the value of the '<em><b>Target Subgraph</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Target Subgraph</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Target Subgraph</em>' containment reference.
   * @see #setTargetSubgraph(Subgraph)
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeTarget_TargetSubgraph()
   * @model containment="true"
   * @generated
   */
  Subgraph getTargetSubgraph();

  /**
   * Sets the value of the '{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget#getTargetSubgraph <em>Target Subgraph</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Target Subgraph</em>' containment reference.
   * @see #getTargetSubgraph()
   * @generated
   */
  void setTargetSubgraph(Subgraph value);

  /**
   * Returns the value of the '<em><b>Targetnode</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Targetnode</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Targetnode</em>' containment reference.
   * @see #setTargetnode(Node)
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeTarget_Targetnode()
   * @model containment="true"
   * @generated
   */
  Node getTargetnode();

  /**
   * Sets the value of the '{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget#getTargetnode <em>Targetnode</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Targetnode</em>' containment reference.
   * @see #getTargetnode()
   * @generated
   */
  void setTargetnode(Node value);

} // EdgeTarget
