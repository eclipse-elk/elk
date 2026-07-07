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
package org.eclipse.elk.core.meta.metaData;

import org.eclipse.emf.common.util.EList;

import org.eclipse.emf.ecore.EObject;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Md Bundle</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdBundle#getLabel <em>Label</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdBundle#getTargetClass <em>Target Class</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdBundle#getDocumentationFolder <em>Documentation Folder</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdBundle#getIdPrefix <em>Id Prefix</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdBundle#getMembers <em>Members</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdBundle()
 * @model
 * @generated
 */
public interface MdBundle extends EObject
{
  /**
   * Returns the value of the '<em><b>Label</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Label</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Label</em>' attribute.
   * @see #setLabel(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdBundle_Label()
   * @model
   * @generated
   */
  String getLabel();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getLabel <em>Label</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Label</em>' attribute.
   * @see #getLabel()
   * @generated
   */
  void setLabel(String value);

  /**
   * Returns the value of the '<em><b>Target Class</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Target Class</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Target Class</em>' attribute.
   * @see #setTargetClass(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdBundle_TargetClass()
   * @model
   * @generated
   */
  String getTargetClass();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getTargetClass <em>Target Class</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Target Class</em>' attribute.
   * @see #getTargetClass()
   * @generated
   */
  void setTargetClass(String value);

  /**
   * Returns the value of the '<em><b>Documentation Folder</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Documentation Folder</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Documentation Folder</em>' attribute.
   * @see #setDocumentationFolder(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdBundle_DocumentationFolder()
   * @model
   * @generated
   */
  String getDocumentationFolder();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getDocumentationFolder <em>Documentation Folder</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Documentation Folder</em>' attribute.
   * @see #getDocumentationFolder()
   * @generated
   */
  void setDocumentationFolder(String value);

  /**
   * Returns the value of the '<em><b>Id Prefix</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Id Prefix</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Id Prefix</em>' attribute.
   * @see #setIdPrefix(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdBundle_IdPrefix()
   * @model
   * @generated
   */
  String getIdPrefix();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdBundle#getIdPrefix <em>Id Prefix</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Id Prefix</em>' attribute.
   * @see #getIdPrefix()
   * @generated
   */
  void setIdPrefix(String value);

  /**
   * Returns the value of the '<em><b>Members</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.core.meta.metaData.MdBundleMember}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Members</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Members</em>' containment reference list.
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdBundle_Members()
   * @model containment="true"
   * @generated
   */
  EList<MdBundleMember> getMembers();

} // MdBundle
