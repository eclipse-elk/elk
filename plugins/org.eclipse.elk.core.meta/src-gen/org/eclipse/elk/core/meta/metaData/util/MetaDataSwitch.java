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
package org.eclipse.elk.core.meta.metaData.util;

import org.eclipse.elk.core.meta.metaData.*;

import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EPackage;

import org.eclipse.emf.ecore.util.Switch;

/**
 * <!-- begin-user-doc -->
 * The <b>Switch</b> for the model's inheritance hierarchy.
 * It supports the call {@link #doSwitch(EObject) doSwitch(object)}
 * to invoke the <code>caseXXX</code> method for each class of the model,
 * starting with the actual class of the object
 * and proceeding up the inheritance hierarchy
 * until a non-null result is returned,
 * which is the result of the switch.
 * <!-- end-user-doc -->
 * @see org.eclipse.elk.core.meta.metaData.MetaDataPackage
 * @generated
 */
public class MetaDataSwitch<T> extends Switch<T>
{
  /**
   * The cached model package
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  protected static MetaDataPackage modelPackage;

  /**
   * Creates an instance of the switch.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @generated
   */
  public MetaDataSwitch()
  {
    if (modelPackage == null)
    {
      modelPackage = MetaDataPackage.eINSTANCE;
    }
  }

  /**
   * Checks whether this is a switch for the given package.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @param ePackage the package in question.
   * @return whether this is a switch for the given package.
   * @generated
   */
  @Override
  protected boolean isSwitchFor(EPackage ePackage)
  {
    return ePackage == modelPackage;
  }

  /**
   * Calls <code>caseXXX</code> for each class of the model until one returns a non null result; it yields that result.
   * <!-- begin-user-doc -->
   * <!-- end-user-doc -->
   * @return the first non-null result returned by a <code>caseXXX</code> call.
   * @generated
   */
  @Override
  protected T doSwitch(int classifierID, EObject theEObject)
  {
    switch (classifierID)
    {
      case MetaDataPackage.MD_MODEL:
      {
        MdModel mdModel = (MdModel)theEObject;
        T result = caseMdModel(mdModel);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_BUNDLE:
      {
        MdBundle mdBundle = (MdBundle)theEObject;
        T result = caseMdBundle(mdBundle);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_BUNDLE_MEMBER:
      {
        MdBundleMember mdBundleMember = (MdBundleMember)theEObject;
        T result = caseMdBundleMember(mdBundleMember);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_GROUP_OR_OPTION:
      {
        MdGroupOrOption mdGroupOrOption = (MdGroupOrOption)theEObject;
        T result = caseMdGroupOrOption(mdGroupOrOption);
        if (result == null) result = caseMdBundleMember(mdGroupOrOption);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_GROUP:
      {
        MdGroup mdGroup = (MdGroup)theEObject;
        T result = caseMdGroup(mdGroup);
        if (result == null) result = caseMdGroupOrOption(mdGroup);
        if (result == null) result = caseMdBundleMember(mdGroup);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_OPTION:
      {
        MdOption mdOption = (MdOption)theEObject;
        T result = caseMdOption(mdOption);
        if (result == null) result = caseMdGroupOrOption(mdOption);
        if (result == null) result = caseMdBundleMember(mdOption);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_OPTION_DEPENDENCY:
      {
        MdOptionDependency mdOptionDependency = (MdOptionDependency)theEObject;
        T result = caseMdOptionDependency(mdOptionDependency);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_ALGORITHM:
      {
        MdAlgorithm mdAlgorithm = (MdAlgorithm)theEObject;
        T result = caseMdAlgorithm(mdAlgorithm);
        if (result == null) result = caseMdBundleMember(mdAlgorithm);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_CATEGORY:
      {
        MdCategory mdCategory = (MdCategory)theEObject;
        T result = caseMdCategory(mdCategory);
        if (result == null) result = caseMdBundleMember(mdCategory);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      case MetaDataPackage.MD_OPTION_SUPPORT:
      {
        MdOptionSupport mdOptionSupport = (MdOptionSupport)theEObject;
        T result = caseMdOptionSupport(mdOptionSupport);
        if (result == null) result = defaultCase(theEObject);
        return result;
      }
      default: return defaultCase(theEObject);
    }
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Model</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Model</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdModel(MdModel object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Bundle</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Bundle</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdBundle(MdBundle object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Bundle Member</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Bundle Member</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdBundleMember(MdBundleMember object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Group Or Option</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Group Or Option</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdGroupOrOption(MdGroupOrOption object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Group</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Group</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdGroup(MdGroup object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Option</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Option</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdOption(MdOption object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Option Dependency</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Option Dependency</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdOptionDependency(MdOptionDependency object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Algorithm</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Algorithm</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdAlgorithm(MdAlgorithm object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Category</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Category</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdCategory(MdCategory object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>Md Option Support</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>Md Option Support</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject) doSwitch(EObject)
   * @generated
   */
  public T caseMdOptionSupport(MdOptionSupport object)
  {
    return null;
  }

  /**
   * Returns the result of interpreting the object as an instance of '<em>EObject</em>'.
   * <!-- begin-user-doc -->
   * This implementation returns null;
   * returning a non-null result will terminate the switch, but this is the last case anyway.
   * <!-- end-user-doc -->
   * @param object the target of the switch.
   * @return the result of interpreting the object as an instance of '<em>EObject</em>'.
   * @see #doSwitch(org.eclipse.emf.ecore.EObject)
   * @generated
   */
  @Override
  public T defaultCase(EObject object)
  {
    return null;
  }

} //MetaDataSwitch
