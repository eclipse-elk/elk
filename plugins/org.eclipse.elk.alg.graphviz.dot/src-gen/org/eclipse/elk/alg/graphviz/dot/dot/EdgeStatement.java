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

import org.eclipse.emf.common.util.EList;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Edge Statement</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement#getSourceNode <em>Source Node</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement#getEdgeTargets <em>Edge Targets</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement#getAttributes <em>Attributes</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeStatement()
 * @model
 * @generated
 */
public interface EdgeStatement extends Statement
{
  /**
   * Returns the value of the '<em><b>Source Node</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Source Node</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Source Node</em>' containment reference.
   * @see #setSourceNode(Node)
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeStatement_SourceNode()
   * @model containment="true"
   * @generated
   */
  Node getSourceNode();

  /**
   * Sets the value of the '{@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeStatement#getSourceNode <em>Source Node</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Source Node</em>' containment reference.
   * @see #getSourceNode()
   * @generated
   */
  void setSourceNode(Node value);

  /**
   * Returns the value of the '<em><b>Edge Targets</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.alg.graphviz.dot.dot.EdgeTarget}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Edge Targets</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Edge Targets</em>' containment reference list.
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeStatement_EdgeTargets()
   * @model containment="true"
   * @generated
   */
  EList<EdgeTarget> getEdgeTargets();

  /**
   * Returns the value of the '<em><b>Attributes</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.alg.graphviz.dot.dot.Attribute}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Attributes</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Attributes</em>' containment reference list.
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getEdgeStatement_Attributes()
   * @model containment="true"
   * @generated
   */
  EList<Attribute> getAttributes();

} // EdgeStatement
