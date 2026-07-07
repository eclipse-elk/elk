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

import org.eclipse.emf.ecore.EFactory;

/**
 * <!-- begin-user-doc -->
 * The <b>Factory</b> for the model.
 * It provides a create method for each non-abstract class of the model.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage
 * @generated
 */
public interface MetaDataFactory extends EFactory
{
  /**
   * The singleton instance of the factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  MetaDataFactory eINSTANCE = org.eclipse.elk.core.meta.metaData.impl.MetaDataFactoryImpl.init();

  /**
   * Returns a new object of class '<em>Md Model</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Model</em>'.
   * @generated
   */
  MdModel createMdModel();

  /**
   * Returns a new object of class '<em>Md Bundle</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Bundle</em>'.
   * @generated
   */
  MdBundle createMdBundle();

  /**
   * Returns a new object of class '<em>Md Bundle Member</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Bundle Member</em>'.
   * @generated
   */
  MdBundleMember createMdBundleMember();

  /**
   * Returns a new object of class '<em>Md Group Or Option</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Group Or Option</em>'.
   * @generated
   */
  MdGroupOrOption createMdGroupOrOption();

  /**
   * Returns a new object of class '<em>Md Group</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Group</em>'.
   * @generated
   */
  MdGroup createMdGroup();

  /**
   * Returns a new object of class '<em>Md Option</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Option</em>'.
   * @generated
   */
  MdOption createMdOption();

  /**
   * Returns a new object of class '<em>Md Option Dependency</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Option Dependency</em>'.
   * @generated
   */
  MdOptionDependency createMdOptionDependency();

  /**
   * Returns a new object of class '<em>Md Algorithm</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Algorithm</em>'.
   * @generated
   */
  MdAlgorithm createMdAlgorithm();

  /**
   * Returns a new object of class '<em>Md Category</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Category</em>'.
   * @generated
   */
  MdCategory createMdCategory();

  /**
   * Returns a new object of class '<em>Md Option Support</em>'.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return a new object of class '<em>Md Option Support</em>'.
   * @generated
   */
  MdOptionSupport createMdOptionSupport();

  /**
   * Returns the package supported by this factory.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the package supported by this factory.
   * @generated
   */
  MetaDataPackage getMetaDataPackage();

} //MetaDataFactory
