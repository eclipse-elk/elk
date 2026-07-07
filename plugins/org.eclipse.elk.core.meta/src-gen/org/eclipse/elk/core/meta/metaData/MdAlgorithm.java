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

import org.eclipse.xtext.common.types.JvmTypeReference;

/**
 * <!-- begin-user-doc -->
 * A representation of the model object '<em><b>Md Algorithm</b></em>'.
 * <!-- end-user-doc -->
 *
 * <p>
 * The following features are supported:
 * </p>
 * <ul>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#isDeprecated <em>Deprecated</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getProvider <em>Provider</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getParameter <em>Parameter</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getLabel <em>Label</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getTargetClass <em>Target Class</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getDescription <em>Description</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getCategory <em>Category</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getPreviewImage <em>Preview Image</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getSupportedFeatures <em>Supported Features</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getValidator <em>Validator</em>}</li>
 *   <li>{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getSupportedOptions <em>Supported Options</em>}</li>
 * </ul>
 *
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm()
 * @model
 * @generated
 */
public interface MdAlgorithm extends MdBundleMember
{
  /**
   * Returns the value of the '<em><b>Deprecated</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Deprecated</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Deprecated</em>' attribute.
   * @see #setDeprecated(boolean)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Deprecated()
   * @model
   * @generated
   */
  boolean isDeprecated();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#isDeprecated <em>Deprecated</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Deprecated</em>' attribute.
   * @see #isDeprecated()
   * @generated
   */
  void setDeprecated(boolean value);

  /**
   * Returns the value of the '<em><b>Provider</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Provider</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Provider</em>' containment reference.
   * @see #setProvider(JvmTypeReference)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Provider()
   * @model containment="true"
   * @generated
   */
  JvmTypeReference getProvider();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getProvider <em>Provider</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Provider</em>' containment reference.
   * @see #getProvider()
   * @generated
   */
  void setProvider(JvmTypeReference value);

  /**
   * Returns the value of the '<em><b>Parameter</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Parameter</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Parameter</em>' attribute.
   * @see #setParameter(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Parameter()
   * @model
   * @generated
   */
  String getParameter();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getParameter <em>Parameter</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Parameter</em>' attribute.
   * @see #getParameter()
   * @generated
   */
  void setParameter(String value);

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
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Label()
   * @model
   * @generated
   */
  String getLabel();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getLabel <em>Label</em>}' attribute.
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
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_TargetClass()
   * @model
   * @generated
   */
  String getTargetClass();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getTargetClass <em>Target Class</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Target Class</em>' attribute.
   * @see #getTargetClass()
   * @generated
   */
  void setTargetClass(String value);

  /**
   * Returns the value of the '<em><b>Description</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Description</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Description</em>' attribute.
   * @see #setDescription(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Description()
   * @model
   * @generated
   */
  String getDescription();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getDescription <em>Description</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Description</em>' attribute.
   * @see #getDescription()
   * @generated
   */
  void setDescription(String value);

  /**
   * Returns the value of the '<em><b>Category</b></em>' reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Category</em>' reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Category</em>' reference.
   * @see #setCategory(MdCategory)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Category()
   * @model
   * @generated
   */
  MdCategory getCategory();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getCategory <em>Category</em>}' reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Category</em>' reference.
   * @see #getCategory()
   * @generated
   */
  void setCategory(MdCategory value);

  /**
   * Returns the value of the '<em><b>Preview Image</b></em>' attribute.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Preview Image</em>' attribute isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Preview Image</em>' attribute.
   * @see #setPreviewImage(String)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_PreviewImage()
   * @model
   * @generated
   */
  String getPreviewImage();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getPreviewImage <em>Preview Image</em>}' attribute.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Preview Image</em>' attribute.
   * @see #getPreviewImage()
   * @generated
   */
  void setPreviewImage(String value);

  /**
   * Returns the value of the '<em><b>Supported Features</b></em>' attribute list.
   * The list contents are of type {@link org.eclipse.elk.core.meta.metaData.MdGraphFeature}.
   * The literals are from the enumeration {@link org.eclipse.elk.core.meta.metaData.MdGraphFeature}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Supported Features</em>' attribute list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Supported Features</em>' attribute list.
   * @see org.eclipse.elk.core.meta.metaData.MdGraphFeature
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_SupportedFeatures()
   * @model unique="false"
   * @generated
   */
  EList<MdGraphFeature> getSupportedFeatures();

  /**
   * Returns the value of the '<em><b>Validator</b></em>' containment reference.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Validator</em>' containment reference isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Validator</em>' containment reference.
   * @see #setValidator(JvmTypeReference)
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_Validator()
   * @model containment="true"
   * @generated
   */
  JvmTypeReference getValidator();

  /**
   * Sets the value of the '{@link org.eclipse.elk.core.meta.metaData.MdAlgorithm#getValidator <em>Validator</em>}' containment reference.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param value the new value of the '<em>Validator</em>' containment reference.
   * @see #getValidator()
   * @generated
   */
  void setValidator(JvmTypeReference value);

  /**
   * Returns the value of the '<em><b>Supported Options</b></em>' containment reference list.
   * The list contents are of type {@link org.eclipse.elk.core.meta.metaData.MdOptionSupport}.
   * <!-- begin-user-doc -->
   * <p>
   * If the meaning of the '<em>Supported Options</em>' containment reference list isn't clear,
   * there really should be more of a description here...
   * </p>
   * <!-- end-user-doc -->
   * @return the value of the '<em>Supported Options</em>' containment reference list.
   * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage#getMdAlgorithm_SupportedOptions()
   * @model containment="true"
   * @generated
   */
  EList<MdOptionSupport> getSupportedOptions();

} // MdAlgorithm
