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
package org.eclipse.elk.core.meta.metaData;

import org.eclipse.emf.ecore.EObject;

import org.eclipse.xtext.xbase.XExpression;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Md Option Support</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getOption <em>Option</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getValue <em>Value</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getDocumentation <em>Documentation</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOptionSupport()
 * @model
 * @generated
 */
public interface MdOptionSupport extends EObject
{
  /**
   * Returns the value of the '<em><b>Option</b></em>' reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Option</em>' reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Option</em>' reference.
   * @see #setOption(MdOption)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOptionSupport_Option()
   * @model
   * @generated
   */
  MdOption getOption();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getOption <em>Option</em>}' reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Option</em>' reference.
   * @see #getOption()
   * @generated
   */
  void setOption(MdOption value);

  /**
   * Returns the value of the '<em><b>Value</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Value</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Value</em>' containment reference.
   * @see #setValue(XExpression)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOptionSupport_Value()
   * @model containment="true"
   * @generated
   */
  XExpression getValue();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getValue <em>Value</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Value</em>' containment reference.
   * @see #getValue()
   * @generated
   */
  void setValue(XExpression value);

  /**
   * Returns the value of the '<em><b>Documentation</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Documentation</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Documentation</em>' attribute.
   * @see #setDocumentation(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdOptionSupport_Documentation()
   * @model
   * @generated
   */
  String getDocumentation();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdOptionSupport#getDocumentation <em>Documentation</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Documentation</em>' attribute.
   * @see #getDocumentation()
   * @generated
   */
  void setDocumentation(String value);

} // MdOptionSupport
