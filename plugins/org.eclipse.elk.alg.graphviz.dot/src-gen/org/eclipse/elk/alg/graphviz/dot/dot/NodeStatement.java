/**
 * ******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.graphviz.dot.dot;

import org.eclipse.emf.common.util.EList;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Node Statement</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement#getNode <em>Node</em>}</li>
 *   <li>{@link org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement#getAttributes <em>Attributes</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getNodeStatement()
 * @model
 * @generated
 */
public interface NodeStatement extends Statement
{
  /**
   * Returns the value of the '<em><b>Node</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Node</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Node</em>' containment reference.
   * @see #setNode(Node)
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getNodeStatement_Node()
   * @model containment="true"
   * @generated
   */
  Node getNode();

  /**
   * Sets the value of the '{@link org.eclipse.elk.alg.graphviz.dot.dot.NodeStatement#getNode <em>Node</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Node</em>' containment reference.
   * @see #getNode()
   * @generated
   */
  void setNode(Node value);

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
   * @see org.eclipse.elk.alg.graphviz.dot.dot.DotPackage#getNodeStatement_Attributes()
   * @model containment="true"
   * @generated
   */
  EList<Attribute> getAttributes();

} // NodeStatement
